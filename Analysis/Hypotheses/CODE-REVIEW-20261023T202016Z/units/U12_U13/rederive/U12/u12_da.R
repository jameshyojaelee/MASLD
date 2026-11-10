#!/usr/bin/env Rscript
# U12 independent re-derivation of the snATAC differential-accessibility counts
# printed in METHODS (39,914 jointly testable; 38,066 indeterminate; 1,844
# source-dependent; 1,505 same sign in the other cohort).
# Written from the METHODS text: donor-by-lineage pseudobulks, MASH vs normal
# donors with >= 20 nuclei, edgeR::filterByExpr with intercept-only design per
# cohort, peaks kept when testable in both cohorts, TMM, voomWithQualityWeights,
# lmFit, robust eBayes, design ~ group, BH within lineage x cohort.
# A peak is supported when q < 0.05 in both cohorts with the same sign.
# Also: (1) leave-one-donor-out sign sensitivity with TMM on the full jointly
# testable family (the producer recomputes TMM on the significant peaks only);
# (2) comparison against the August r1 and 2026-10-22 corrected tables.

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(edgeR); library(limma)
})
set.seed(1)  # nothing below draws random numbers; set for the record

args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CAND <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates"
R1 <- file.path(CAND, "atac-context-v3-candidate-2026-08-11-r1")
CORR <- file.path(CAND, "atac-context-v3-statistics-correction-20261022")
COHORTS <- c("GSE244832", "GSE281367")
LINEAGES <- c("hepatocyte", "stellate", "macrophage")

read_lineage <- function(cohort, lineage) {
  base <- file.path(R1, "counts", cohort)
  con <- gzfile(file.path(base, paste0(lineage, ".mtx.gz")), "rb")
  m <- as(readMM(con), "CsparseMatrix"); close(con)
  donors <- fread(file.path(base, paste0(lineage, ".donors.tsv")))
  peaks <- fread(file.path(base, paste0(lineage, ".peaks.tsv")))
  stopifnot(nrow(m) == nrow(donors), ncol(m) == nrow(peaks))
  dimnames(m) <- list(donors$donor_id, peaks$peak_coordinate)
  list(m = m, donors = donors)
}

# Eligibility from METHODS: MASH or normal, >= 20 nuclei in the lineage.
eligible_donors <- function(d) d[condition %in% c("NORMAL", "MASH") & n_cells >= 20]

fit_cohort <- function(counts_dp, donors) {
  # counts_dp: donors x peaks
  grp <- factor(donors$condition, levels = c("NORMAL", "MASH"))
  X <- model.matrix(~ grp)
  dge <- DGEList(counts = t(as.matrix(counts_dp)))
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voomWithQualityWeights(dge, X)
  fit <- eBayes(lmFit(v, X), robust = TRUE)
  p <- fit$p.value[, 2]
  data.table(
    peak_coordinate = rownames(fit$coefficients),
    logFC = fit$coefficients[, 2],
    se_moderated = fit$stdev.unscaled[, 2] * sqrt(fit$s2.post),
    se_sigma = fit$stdev.unscaled[, 2] * fit$sigma,
    t = fit$t[, 2],
    p = p,
    q = p.adjust(p, "BH"),
    df_total = fit$df.total
  )
}

lodo_signs <- function(counts_dp, donors, peaks_eval, ref_logfc) {
  # TMM and voom on the full jointly testable family; signs read on peaks_eval.
  frac <- matrix(NA_real_, length(peaks_eval), nrow(donors))
  for (i in seq_len(nrow(donors))) {
    keep <- seq_len(nrow(donors)) != i
    tab <- table(donors$condition[keep])
    if (any(tab[c("NORMAL", "MASH")] < 3)) next
    grp <- factor(donors$condition[keep], levels = c("NORMAL", "MASH"))
    X <- model.matrix(~ grp)
    dge <- calcNormFactors(DGEList(counts = t(as.matrix(counts_dp[keep, , drop = FALSE]))), method = "TMM")
    v <- voom(dge, X)
    cf <- lmFit(v, X)$coefficients[peaks_eval, 2]
    frac[, i] <- sign(cf) == sign(ref_logfc)
  }
  rowMeans(frac, na.rm = TRUE)
}

donor_rows <- list(); summary_rows <- list(); joined_all <- list()
for (lin in LINEAGES) {
  x <- lapply(COHORTS, read_lineage, lineage = lin); names(x) <- COHORTS
  keep_peaks <- list(); el <- list()
  for (co in COHORTS) {
    d <- eligible_donors(x[[co]]$donors)
    el[[co]] <- d
    donor_rows[[paste(lin, co)]] <- data.table(
      lineage = lin, cohort = co,
      n_normal = sum(d$condition == "NORMAL"), n_mash = sum(d$condition == "MASH"),
      n_rows_matrix = nrow(x[[co]]$donors),
      eligible_matches_contrast_eligible = identical(
        sort(d$donor_id),
        sort(x[[co]]$donors[contrast_eligible == TRUE & condition %in% c("NORMAL", "MASH"), donor_id])),
      unique_donor_ids = !anyDuplicated(x[[co]]$donors$donor_id),
      min_nuclei_eligible = min(d$n_cells))
    cm <- x[[co]]$m[d$donor_id, , drop = FALSE]
    dge0 <- DGEList(counts = t(as.matrix(cm)))
    keep_peaks[[co]] <- colnames(cm)[filterByExpr(dge0, design = matrix(1, ncol(dge0), 1))]
  }
  common <- intersect(keep_peaks[[1]], keep_peaks[[2]])
  if (length(common) == 0) {
    summary_rows[[lin]] <- data.table(lineage = lin, n_filter_244 = length(keep_peaks[[1]]),
      n_filter_281 = length(keep_peaks[[2]]), joint = 0L)
    next
  }
  fits <- lapply(COHORTS, function(co) {
    fit_cohort(x[[co]]$m[el[[co]]$donor_id, common, drop = FALSE], el[[co]])
  })
  names(fits) <- COHORTS
  j <- merge(fits[[1]], fits[[2]], by = "peak_coordinate", suffixes = c("_244", "_281"))
  j[, sig244 := q_244 < 0.05][, sig281 := q_281 < 0.05]
  j[, state := fifelse(sig244 & sig281,
                       fifelse(sign(logFC_244) == sign(logFC_281), "supported", "discordant"),
                       fifelse(sig244 | sig281, "source_dependent", "indeterminate"))]
  j[, same_sign := sign(logFC_244) == sign(logFC_281)]
  # LODO on peaks significant in either cohort, TMM over the full family.
  cand <- j[sig244 | sig281, peak_coordinate]
  j[, lodo_full_244 := NA_real_][, lodo_full_281 := NA_real_]
  if (length(cand)) {
    for (co in COHORTS) {
      ref <- j[match(cand, peak_coordinate), get(paste0("logFC_", substr(co, 4, 6)))]
      fr <- lodo_signs(x[[co]]$m[el[[co]]$donor_id, common, drop = FALSE], el[[co]], cand, ref)
      j[match(cand, peak_coordinate), (paste0("lodo_full_", substr(co, 4, 6))) := fr]
    }
  }
  j[, lineage := lin]
  joined_all[[lin]] <- j
  summary_rows[[lin]] <- data.table(
    lineage = lin, n_filter_244 = length(keep_peaks[[1]]), n_filter_281 = length(keep_peaks[[2]]),
    joint = nrow(j), sig244 = sum(j$sig244), sig281 = sum(j$sig281),
    supported = sum(j$state == "supported"), discordant = sum(j$state == "discordant"),
    source_dependent = sum(j$state == "source_dependent"),
    source_dependent_same_sign = sum(j$state == "source_dependent" & j$same_sign),
    indeterminate = sum(j$state == "indeterminate"))
}
donors_tab <- rbindlist(donor_rows); summ <- rbindlist(summary_rows, fill = TRUE)
mine <- rbindlist(joined_all)
tot <- summ[, lapply(.SD, sum, na.rm = TRUE), .SDcols = setdiff(names(summ), "lineage")][, lineage := "TOTAL"]
summ <- rbind(summ, tot, fill = TRUE)
fwrite(donors_tab, file.path(OUT, "da_donors.tsv"), sep = "\t")
fwrite(summ, file.path(OUT, "da_state_counts.tsv"), sep = "\t")
fwrite(mine, file.path(OUT, "da_peak_rederived.tsv.gz"), sep = "\t")

# Compare with producer tables.
r1 <- fread(file.path(R1, "da/da_peak_results.tsv.gz"))
co <- fread(file.path(CORR, "da/da_peak_results.tsv.gz"))
cmpr <- function(prod, tag) {
  m <- merge(mine, prod, by = c("lineage", "peak_coordinate"), all = TRUE)
  data.table(
    table = tag,
    rows_mine = nrow(mine), rows_prod = nrow(prod), rows_joined = nrow(m),
    only_mine = m[is.na(evidence_state), .N], only_prod = m[is.na(state), .N],
    max_abs_logFC_244 = m[, max(abs(logFC_244 - logFC_gse244832), na.rm = TRUE)],
    max_abs_logFC_281 = m[, max(abs(logFC_281 - logFC_gse281367), na.rm = TRUE)],
    max_abs_p_244 = m[, max(abs(p_244 - pvalue_gse244832), na.rm = TRUE)],
    max_abs_p_281 = m[, max(abs(p_281 - pvalue_gse281367), na.rm = TRUE)],
    max_abs_q_244 = m[, max(abs(q_244 - qvalue_gse244832), na.rm = TRUE)],
    max_abs_q_281 = m[, max(abs(q_281 - qvalue_gse281367), na.rm = TRUE)],
    max_abs_se_moderated_244 = m[, max(abs(se_moderated_244 - standard_error_gse244832), na.rm = TRUE)],
    max_abs_se_moderated_281 = m[, max(abs(se_moderated_281 - standard_error_gse281367), na.rm = TRUE)],
    max_abs_se_sigma_244 = m[, max(abs(se_sigma_244 - standard_error_gse244832), na.rm = TRUE)],
    max_abs_se_sigma_281 = m[, max(abs(se_sigma_281 - standard_error_gse281367), na.rm = TRUE)],
    state_disagreements = m[!is.na(state) & !is.na(evidence_state) & state != evidence_state, .N],
    producer_lodo_flag_false_among_source_dependent = m[evidence_state == "source_dependent" & sensitivity_no_reversal == FALSE, .N],
    producer_lodo_flag_false_among_supported = m[evidence_state == "supported" & sensitivity_no_reversal == FALSE, .N]
  )
}
cmp <- rbind(cmpr(r1, "r1_2026-08-11"), cmpr(co, "corrected_2026-10-22"))
fwrite(cmp, file.path(OUT, "da_compare_producer.tsv"), sep = "\t")

# Sensitivity rule with LODO computed on full-family TMM (>= 0.80 same sign),
# applied to supported and source-dependent peaks.
m <- merge(mine, r1[, .(lineage, peak_coordinate, lodo_sign_fraction_gse244832, lodo_sign_fraction_gse281367,
                        batch_adjusted_logFC_gse244832, sensitivity_no_reversal)],
           by = c("lineage", "peak_coordinate"))
m[, my_stable := (is.na(lodo_full_244) | lodo_full_244 >= 0.8) & (is.na(lodo_full_281) | lodo_full_281 >= 0.8) &
    (is.na(batch_adjusted_logFC_gse244832) | sign(batch_adjusted_logFC_gse244832) == sign(logFC_244))]
sens <- m[state %in% c("supported", "source_dependent"), .(
  n = .N,
  producer_flag_false = sum(!sensitivity_no_reversal),
  fulltmm_flag_false = sum(!my_stable),
  flags_disagree = sum(my_stable != sensitivity_no_reversal),
  max_abs_lodo_diff_244 = max(abs(lodo_full_244 - lodo_sign_fraction_gse244832), na.rm = TRUE),
  max_abs_lodo_diff_281 = max(abs(lodo_full_281 - lodo_sign_fraction_gse281367), na.rm = TRUE)
), by = .(lineage, state)]
fwrite(sens, file.path(OUT, "da_sensitivity_rule.tsv"), sep = "\t")

print(donors_tab); print(summ); print(cmp); print(sens)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
