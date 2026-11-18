#!/usr/bin/env Rscript
# Donor-level differential accessibility (20261118 candidate, N10/N19).
#
# Key methodological updates over 04_run_da.R:
# 1. Exclude both conflicted donors D11 (MASH) and D14 (MASL) from GSE244832 before
#    fitting (yielding 5 NORMAL vs 8 MASH in primary GSE244832 contrast), and record
#    the explicit exclusion audit in atac_da_conflict_audit_20261118.tsv.
# 2. Fit `batch_mm`-adjusted GSE244832 voomWithQualityWeights (`~ factor(batch_mm) + group`)
#    as primary/co-primary contrast alongside unadjusted (`~ group`).
# 3. Replace `binom.test(n_same, n_discovery, p=0.5)` with Freedman-Lane empirical
#    donor-residual permutation p-values (N_FL = 9,999 draws) and exact within-batch_mm
#    donor sign-flip p-values (2^5 = 32 balanced batch splits), with BH q-values across
#    testable lineages (`directional_empirical_qvalue`), respecting inter-peak correlation
#    and biological donor units.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(GenomicRanges)
  library(rtracklayer)
  library(parallel)
})

set.seed(20261118)
N_CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
N_FL <- 9999L

BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RELEASE_ID <- "atac-context-v3-candidate-2026-08-11-r1"
CANDIDATE <- file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID)
OUT <- Sys.getenv(
  "ATAC_DA_OUT_20261118",
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates/atac_da_20261118")
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GTF <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
G244_META <- file.path(BASE, "data/GSE244832/metadata/donor_pairing.csv")
LINEAGES <- c("hepatocyte", "stellate", "macrophage")
CONFLICT_DONORS <- c("D11", "D14")

read_sparse_gz <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  as(readMM(con), "dgCMatrix")
}

load_counts <- function(cohort, lineage, root_override = NULL) {
  root <- if (is.null(root_override)) file.path(CANDIDATE, "counts", cohort) else file.path(root_override, cohort)
  counts <- read_sparse_gz(file.path(root, paste0(lineage, ".mtx.gz")))
  donors <- fread(file.path(root, paste0(lineage, ".donors.tsv")))
  peaks <- fread(file.path(root, paste0(lineage, ".peaks.tsv")))
  if (nrow(counts) != nrow(donors) || ncol(counts) != nrow(peaks)) stop("Sparse count dimension drift for ", cohort, "/", lineage)
  rownames(counts) <- donors$donor_id
  colnames(counts) <- peaks$peak_coordinate
  list(counts = counts, donors = donors, peaks = peaks)
}

eligible_primary <- function(x, exclude_ids = character()) {
  keep <- x$donors$contrast_eligible == TRUE &
    x$donors$condition %chin% c("NORMAL", "MASH") &
    !(x$donors$donor_id %chin% exclude_ids)
  counts <- x$counts[keep, , drop = FALSE]
  donors <- x$donors[keep]
  observed <- donors[, .N, by = condition]
  n_normal <- observed[condition == "NORMAL", N]; if (length(n_normal) == 0L) n_normal <- 0L
  n_mash <- observed[condition == "MASH", N]; if (length(n_mash) == 0L) n_mash <- 0L
  list(counts = counts, donors = donors, n_normal = n_normal, n_mash = n_mash)
}

condition_blind_filter <- function(counts) {
  y <- DGEList(counts = t(counts))
  design <- matrix(1, nrow = ncol(y), ncol = 1L)
  filterByExpr(y, design = design)
}

fit_voom_full <- function(counts, donors, batch = FALSE) {
  group <- factor(donors$condition, levels = c("NORMAL", "MASH"))
  if (batch) {
    batch_factor <- factor(donors$batch_mm)
    design <- model.matrix(~ batch_factor + group)
    coefficient <- ncol(design)
  } else {
    design <- model.matrix(~ group)
    coefficient <- 2L
  }
  if (qr(design)$rank != ncol(design)) stop("DA design is not full rank")
  y <- DGEList(counts = t(counts))
  y <- calcNormFactors(y, method = "TMM")
  v <- voomWithQualityWeights(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  table <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  answer <- data.table(
    peak_coordinate = rownames(table),
    logFC = table$logFC,
    AveExpr = table$AveExpr,
    t = table$t,
    standard_error = sqrt(fit$s2.post) * fit$stdev.unscaled[, coefficient],
    pvalue = table$P.Value,
    qvalue = p.adjust(table$P.Value, method = "BH"),
    B = table$B
  )
  qc <- data.table(
    donor_id = donors$donor_id,
    condition = as.character(group),
    raw_library_size = y$samples$lib.size,
    tmm_normalization_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors,
    sample_quality_weight = v$targets$sample.weights,
    design_rank = qr(design)$rank,
    n_model_columns = ncol(design),
    n_tested_peaks = nrow(answer)
  )
  list(table = answer, qc = qc, v = v, design = design, coefficient = coefficient, donors = donors)
}

fit_simple_effect <- function(counts, donors, batch = FALSE) {
  group <- factor(donors$condition, levels = c("NORMAL", "MASH"))
  if (batch) {
    batch_factor <- factor(donors$batch_mm)
    design <- model.matrix(~ batch_factor + group)
    coef_idx <- ncol(design)
  } else {
    design <- model.matrix(~ group)
    coef_idx <- 2L
  }
  if (qr(design)$rank != ncol(design)) stop("Sensitivity design is not full rank")
  y <- calcNormFactors(DGEList(counts = t(counts)), method = "TMM")
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  setNames(as.numeric(fit$coefficients[, coef_idx]), rownames(fit$coefficients))
}

# Freedman-Lane donor-relabeling test holding GSE281367 fixed and permuting
# reduced-model (~ batch_mm) residuals of GSE244832 (or reduced intercept for unadjusted)
freedman_lane_da <- function(fit_obj, disc_idx, sign_281_disc, obs_logFC_281, n_same_obs, r_obs, seed) {
  E <- fit_obj$v$E           # P x N
  W <- fit_obj$v$weights     # P x N
  X <- fit_obj$design        # N x p (last column is groupMASH)
  Z <- X[, -fit_obj$coefficient, drop = FALSE]
  P <- nrow(E); N <- ncol(E)
  A <- matrix(0, nrow = P, ncol = N)
  E0 <- matrix(0, nrow = P, ncol = N)
  for (g in seq_len(P)) {
    wg <- W[g, ]
    XtWX_inv <- solve(crossprod(X, wg * X))
    A[g, ] <- (XtWX_inv[fit_obj$coefficient, ] %*% t(X)) * wg
    ZtWZ_inv <- solve(crossprod(Z, wg * Z))
    gamma0 <- ZtWZ_inv %*% crossprod(Z, wg * E[g, ])
    E0[g, ] <- E[g, ] - as.numeric(Z %*% gamma0)
  }
  max_diff <- max(abs(rowSums(A * E0) - fit_obj$table$logFC))
  if (max_diff > 1e-8) stop("Freedman-Lane linear projection mismatch: ", max_diff)

  set.seed(seed)
  perms <- lapply(seq_len(N_FL), function(b) sample.int(N))
  chunks <- split(seq_len(N_FL), cut(seq_len(N_FL), breaks = min(N_CORES, N_FL), labels = FALSE))
  res_list <- mclapply(chunks, function(idxs) {
    ns_vec <- integer(length(idxs))
    rc_vec <- numeric(length(idxs))
    for (k in seq_along(idxs)) {
      p <- perms[[idxs[k]]]
      b_star <- rowSums(A * E0[, p, drop = FALSE])
      ns_vec[k] <- if (length(disc_idx) > 0L) sum(sign(b_star[disc_idx]) == sign_281_disc) else 0L
      rc_vec[k] <- cor(b_star, obs_logFC_281, method = "pearson")
    }
    list(ns = ns_vec, rc = rc_vec)
  }, mc.cores = min(N_CORES, length(chunks)))
  ns_all <- unlist(lapply(res_list, `[[`, "ns"))
  rc_all <- unlist(lapply(res_list, `[[`, "rc"))
  p_fl_conc <- if (length(disc_idx) > 0L) (1 + sum(ns_all >= n_same_obs)) / (N_FL + 1) else 1.0
  p_fl_corr <- (1 + sum(rc_all >= r_obs - 1e-12, na.rm = TRUE)) / (N_FL + 1)

  # Also exact within-batch_mm balanced sign-flip (2^5 = 32 splits for batches with 1N+1M)
  d_tab <- fit_obj$donors
  paired_batches <- names(which(tapply(d_tab$condition, d_tab$batch_mm, function(cc) sum(cc == "NORMAL") == 1 && sum(cc == "MASH") == 1)))
  p_wb <- NA_real_
  if (length(paired_batches) >= 1L) {
    pair_idx <- lapply(paired_batches, function(b) which(as.character(d_tab$batch_mm) == b))
    grid_flips <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(pair_idx))))
    ns_wb <- integer(nrow(grid_flips))
    for (r_i in seq_len(nrow(grid_flips))) {
      p <- seq_len(N)
      for (b_i in seq_along(pair_idx)) {
        if (grid_flips[r_i, b_i]) p[pair_idx[[b_i]]] <- rev(pair_idx[[b_i]])
      }
      b_star <- rowSums(A * E0[, p, drop = FALSE])
      ns_wb[r_i] <- if (length(disc_idx) > 0L) sum(sign(b_star[disc_idx]) == sign_281_disc) else 0L
    }
    p_wb <- if (length(disc_idx) > 0L) mean(ns_wb >= n_same_obs) else 1.0
  }
  list(p_fl_conc = p_fl_conc, p_fl_corr = p_fl_corr, p_within_batch_exact = p_wb, n_wb_splits = ifelse(length(paired_batches) > 0, 2^length(paired_batches), 0L))
}

parse_peaks <- function(values) {
  match <- regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", values)
  pieces <- regmatches(values, match)
  data.table(
    peak_coordinate = values,
    chrom = vapply(pieces, `[`, character(1), 2L),
    start0 = as.integer(vapply(pieces, `[`, character(1), 3L)),
    end = as.integer(vapply(pieces, `[`, character(1), 4L))
  )
}

message("[DA_20261118] Importing GENCODE v49 promoters")
genes <- import(GTF, format = "gtf", feature.type = "gene")
genes <- genes[as.character(seqnames(genes)) %chin% paste0("chr", c(1:22, "X", "Y"))]
names(genes) <- as.character(mcols(genes)$gene_name)
genes <- genes[!is.na(names(genes)) & names(genes) != ""]
promoters_v49 <- promoters(genes, upstream = 2000L, downstream = 2001L)

annotate_promoters <- function(peaks) {
  parsed <- parse_peaks(peaks)
  ranges <- GRanges(seqnames = parsed$chrom, ranges = IRanges(start = parsed$start0 + 1L, end = parsed$end))
  hits <- findOverlaps(ranges, promoters_v49, ignore.strand = TRUE)
  by_peak <- split(names(promoters_v49)[subjectHits(hits)], queryHits(hits))
  annotation <- rep("", length(peaks))
  for (index in names(by_peak)) annotation[as.integer(index)] <- paste(sort(unique(by_peak[[index]])), collapse = ";")
  annotation
}

g244_metadata <- fread(G244_META, select = c("donor_id", "batch_mm", "condition", "confidence"))
fwrite(g244_metadata[, .(donor_id, condition, batch_mm, confidence, excluded_20261118 = donor_id %chin% CONFLICT_DONORS)],
       file.path(OUT, "atac_da_conflict_audit_20261118.tsv"), sep = "\t")

joined_results <- list()
lineage_summary <- list()
model_qc_rows <- list()
filter_qc_rows <- list()

for (l_idx in seq_along(LINEAGES)) {
  lineage <- LINEAGES[l_idx]
  message("[DA_20261118] ", lineage)
  raw244 <- load_counts("GSE244832", lineage)
  raw281 <- load_counts("GSE281367", lineage)
  raw244$donors <- merge(raw244$donors, g244_metadata[, .(donor_id, batch_mm)], by = "donor_id", all.x = TRUE, sort = FALSE)
  raw244$donors <- raw244$donors[match(rownames(raw244$counts), donor_id)]

  x244_clean <- eligible_primary(raw244, exclude_ids = CONFLICT_DONORS)
  x244_orig <- eligible_primary(raw244, exclude_ids = character())
  x281 <- eligible_primary(raw281)

  keep244 <- condition_blind_filter(x244_clean$counts)
  keep281 <- condition_blind_filter(x281$counts)
  common <- intersect(colnames(x244_clean$counts)[keep244], colnames(x281$counts)[keep281])

  filter_qc_rows[[lineage]] <- data.table(
    lineage = lineage,
    n_normal_gse244832_clean = x244_clean$n_normal,
    n_mash_gse244832_clean = x244_clean$n_mash,
    n_normal_gse244832_orig = x244_orig$n_normal,
    n_mash_gse244832_orig = x244_orig$n_mash,
    excluded_conflicts = paste(CONFLICT_DONORS, collapse = ","),
    n_filter_retained_gse244832 = sum(keep244),
    n_filter_retained_gse281367 = sum(keep281),
    n_jointly_testable = length(common)
  )

  if (length(common) == 0L) {
    lineage_summary[[lineage]] <- data.table(
      lineage = lineage, lineage_state = "untestable",
      untestable_reason = "no_jointly_testable_peaks_after_filterByExpr",
      n_normal_gse244832 = x244_clean$n_normal, n_mash_gse244832 = x244_clean$n_mash,
      n_normal_gse281367 = x281$n_normal, n_mash_gse281367 = x281$n_mash,
      n_jointly_testable = 0L, n_sig_gse244832_batch_adj = 0L, n_sig_gse244832_unadj = 0L,
      n_sig_gse281367 = 0L, n_supported_batch_adj = 0L, n_supported_unadj = 0L,
      directional_empirical_pvalue = 1.0
    )
    next
  }

  counts244 <- x244_clean$counts[, common, drop = FALSE]
  counts281 <- x281$counts[, common, drop = FALSE]

  fit244_adj_obj <- fit_voom_full(counts244, x244_clean$donors, batch = TRUE)
  fit244_unadj_obj <- fit_voom_full(counts244, x244_clean$donors, batch = FALSE)
  fit281_obj <- fit_voom_full(counts281, x281$donors, batch = FALSE)

  fit244_adj <- fit244_adj_obj$table
  fit244_unadj <- fit244_unadj_obj$table
  fit281 <- fit281_obj$table

  fit244_adj_obj$qc[, `:=`(cohort = "GSE244832", lineage = lineage, model = "primary_batch_mm_adjusted_D11D14_excl")]
  fit244_unadj_obj$qc[, `:=`(cohort = "GSE244832", lineage = lineage, model = "coprimary_unadjusted_D11D14_excl")]
  fit281_obj$qc[, `:=`(cohort = "GSE281367", lineage = lineage, model = "primary_condition_only")]
  model_qc_rows[[paste(lineage, "244_adj")]] <- fit244_adj_obj$qc
  model_qc_rows[[paste(lineage, "244_unadj")]] <- fit244_unadj_obj$qc
  model_qc_rows[[paste(lineage, "281")]] <- fit281_obj$qc

  fwrite(fit244_adj, file.path(OUT, paste0("cohort_primary_batch_adj_GSE244832_", lineage, ".tsv")), sep = "\t")
  fwrite(fit244_unadj, file.path(OUT, paste0("cohort_coprimary_unadj_GSE244832_", lineage, ".tsv")), sep = "\t")
  fwrite(fit281, file.path(OUT, paste0("cohort_primary_GSE281367_", lineage, ".tsv")), sep = "\t")

  joined <- merge(fit244_adj, fit281, by = "peak_coordinate", suffixes = c("_gse244832", "_gse281367"))
  joined <- merge(
    joined,
    fit244_unadj[, .(peak_coordinate, unadj_logFC_gse244832 = logFC, unadj_pvalue_gse244832 = pvalue, unadj_qvalue_gse244832 = qvalue)],
    by = "peak_coordinate"
  )
  joined[, evidence_state := fifelse(
    qvalue_gse244832 < 0.05 & qvalue_gse281367 < 0.05,
    fifelse(logFC_gse244832 * logFC_gse281367 > 0, "supported", "discordant"),
    fifelse(xor(qvalue_gse244832 < 0.05, qvalue_gse281367 < 0.05), "source_dependent", "indeterminate")
  )]
  joined[, unadj_evidence_state := fifelse(
    unadj_qvalue_gse244832 < 0.05 & qvalue_gse281367 < 0.05,
    fifelse(unadj_logFC_gse244832 * logFC_gse281367 > 0, "supported", "discordant"),
    fifelse(xor(unadj_qvalue_gse244832 < 0.05, qvalue_gse281367 < 0.05), "source_dependent", "indeterminate")
  )]
  joined[, promoter_genes := annotate_promoters(peak_coordinate)]

  # Original 5 vs 9 (D11 included, unadjusted) sensitivity logFC for audit comparison
  keep_orig_common <- intersect(common, colnames(x244_orig$counts))
  orig_eff <- fit_simple_effect(x244_orig$counts[, keep_orig_common, drop = FALSE], x244_orig$donors, batch = FALSE)
  joined[, orig_d11incl_unadj_logFC_gse244832 := orig_eff[peak_coordinate]]

  # Broad MASL+MASH excluding D11 & D14
  broad_keep <- raw244$donors$contrast_eligible == TRUE &
    raw244$donors$condition %chin% c("NORMAL", "MASL", "MASH") &
    !(raw244$donors$donor_id %chin% CONFLICT_DONORS)
  broad_donors <- copy(raw244$donors[broad_keep])
  broad_donors[, condition := fifelse(condition == "NORMAL", "NORMAL", "MASH")]
  broad_eff <- fit_simple_effect(raw244$counts[broad_keep, common, drop = FALSE], broad_donors, batch = TRUE)
  joined[, broad_masld_batch_adj_logFC_gse244832 := broad_eff[peak_coordinate]]
  joined[, lineage := lineage]
  joined_results[[lineage]] <- joined

  # Order-aligned vectors for Freedman-Lane
  idx_281 <- match(fit244_adj$peak_coordinate, fit281$peak_coordinate)
  stopifnot(!anyNA(idx_281))
  lfc281_vec <- fit281$logFC[idx_281]
  q281_vec <- fit281$qvalue[idx_281]
  disc_idx <- which(q281_vec < 0.05)
  sign_281_disc <- sign(lfc281_vec[disc_idx])

  n_disc <- length(disc_idx)
  n_same_adj <- if (n_disc > 0L) sum(sign(fit244_adj$logFC[disc_idx]) == sign_281_disc) else 0L
  n_same_unadj <- if (n_disc > 0L) sum(sign(fit244_unadj$logFC[disc_idx]) == sign_281_disc) else 0L
  r_pear_adj <- cor(fit244_adj$logFC, lfc281_vec, method = "pearson")
  r_spear_adj <- cor(fit244_adj$logFC, lfc281_vec, method = "spearman")
  r_pear_unadj <- cor(fit244_unadj$logFC, lfc281_vec, method = "pearson")

  fl_adj <- freedman_lane_da(fit244_adj_obj, disc_idx, sign_281_disc, lfc281_vec, n_same_adj, r_pear_adj, 20261118L + 100L * l_idx)
  fl_unadj <- freedman_lane_da(fit244_unadj_obj, disc_idx, sign_281_disc, lfc281_vec, n_same_unadj, r_pear_unadj, 20262118L + 100L * l_idx)

  lineage_summary[[lineage]] <- data.table(
    lineage = lineage,
    lineage_state = "testable",
    untestable_reason = "",
    n_normal_gse244832 = x244_clean$n_normal,
    n_mash_gse244832 = x244_clean$n_mash,
    n_normal_gse281367 = x281$n_normal,
    n_mash_gse281367 = x281$n_mash,
    n_jointly_testable = nrow(joined),
    n_sig_gse244832_batch_adj = sum(joined$qvalue_gse244832 < 0.05),
    n_sig_gse244832_unadj = sum(joined$unadj_qvalue_gse244832 < 0.05),
    n_sig_gse281367 = sum(joined$qvalue_gse281367 < 0.05),
    n_supported_batch_adj = sum(joined$evidence_state == "supported"),
    n_discordant_batch_adj = sum(joined$evidence_state == "discordant"),
    n_supported_unadj = sum(joined$unadj_evidence_state == "supported"),
    n_discordant_unadj = sum(joined$unadj_evidence_state == "discordant"),
    pearson_effect_correlation_batch_adj = r_pear_adj,
    spearman_effect_correlation_batch_adj = r_spear_adj,
    pearson_effect_correlation_unadj = r_pear_unadj,
    n_gse281367_discovery = n_disc,
    n_direction_concordant_batch_adj = n_same_adj,
    n_direction_concordant_unadj = n_same_unadj,
    directional_empirical_pvalue = fl_adj$p_fl_conc,
    correlation_empirical_pvalue = fl_adj$p_fl_corr,
    within_batch_exact32_pvalue = fl_adj$p_within_batch_exact,
    unadj_directional_empirical_pvalue = fl_unadj$p_fl_conc,
    unadj_correlation_empirical_pvalue = fl_unadj$p_fl_corr
  )
}

all_joined <- rbindlist(joined_results, fill = TRUE)
fwrite(all_joined, file.path(OUT, "da_peak_results.tsv.gz"), sep = "\t")
fwrite(rbindlist(filter_qc_rows), file.path(OUT, "primary_peak_filter_qc.tsv"), sep = "\t")

summary <- rbindlist(lineage_summary, fill = TRUE)
summary[lineage_state == "testable", directional_empirical_qvalue := p.adjust(directional_empirical_pvalue, method = "BH")]
summary[lineage_state == "testable", correlation_empirical_qvalue := p.adjust(correlation_empirical_pvalue, method = "BH")]
summary[lineage_state == "testable", unadj_directional_empirical_qvalue := p.adjust(unadj_directional_empirical_pvalue, method = "BH")]
fwrite(summary, file.path(OUT, "da_lineage_summary.tsv"), sep = "\t")
fwrite(rbindlist(model_qc_rows, fill = TRUE), file.path(OUT, "da_model_qc.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(summary)
message("[DA_20261118] Done -> ", OUT)
