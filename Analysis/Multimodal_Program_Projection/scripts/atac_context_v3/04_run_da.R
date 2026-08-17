#!/usr/bin/env Rscript

# Donor-level, separate-cohort differential accessibility on exact consensus counts.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(GenomicRanges)
  library(rtracklayer)
})

set.seed(42)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "atac-context-v3-candidate-2026-08-11-r1"
CANDIDATE <- Sys.getenv(
  "ATAC_V3_CANDIDATE_ROOT",
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID)
)
EXPECTED <- normalizePath(
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID),
  mustWork = FALSE
)
if (!identical(normalizePath(CANDIDATE, mustWork = FALSE), EXPECTED)) {
  stop("Unsafe candidate root: ", CANDIDATE)
}
if (!file.exists(file.path(CANDIDATE, "consensus_peak_manifest.tsv"))) {
  stop("Consensus peak stage is incomplete")
}
OUT <- file.path(CANDIDATE, "da")
if (dir.exists(OUT)) stop("Refusing to overwrite DA output: ", OUT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GTF <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
G244_META <- file.path(BASE, "data/GSE244832/metadata/donor_pairing.csv")
LINEAGES <- c("hepatocyte", "stellate", "macrophage")

read_sparse_gz <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  as(readMM(con), "dgCMatrix")
}

load_counts <- function(cohort, lineage, root_override = NULL) {
  root <- if (is.null(root_override)) {
    file.path(CANDIDATE, "counts", cohort)
  } else {
    file.path(root_override, cohort)
  }
  matrix_path <- file.path(root, paste0(lineage, ".mtx.gz"))
  donors_path <- file.path(root, paste0(lineage, ".donors.tsv"))
  peaks_path <- file.path(root, paste0(lineage, ".peaks.tsv"))
  if (!all(file.exists(matrix_path, donors_path, peaks_path))) {
    stop("Missing exact count artifact for ", cohort, "/", lineage)
  }
  counts <- read_sparse_gz(matrix_path)
  donors <- fread(donors_path)
  peaks <- fread(peaks_path)
  if (nrow(counts) != nrow(donors) || ncol(counts) != nrow(peaks)) {
    stop("Sparse count dimension drift for ", cohort, "/", lineage)
  }
  if (anyDuplicated(peaks$peak_coordinate)) stop("Duplicated peak coordinates")
  rownames(counts) <- donors$donor_id
  colnames(counts) <- peaks$peak_coordinate
  list(counts = counts, donors = donors, peaks = peaks)
}

eligible_primary <- function(x) {
  keep <- x$donors$contrast_eligible == TRUE & x$donors$condition %chin% c("NORMAL", "MASH")
  counts <- x$counts[keep, , drop = FALSE]
  donors <- x$donors[keep]
  observed <- donors[, .N, by = condition]
  n_normal <- observed[condition == "NORMAL", N]
  n_mash <- observed[condition == "MASH", N]
  if (length(n_normal) == 0L) n_normal <- 0L
  if (length(n_mash) == 0L) n_mash <- 0L
  list(counts = counts, donors = donors, n_normal = n_normal, n_mash = n_mash)
}

condition_blind_filter <- function(counts) {
  y <- DGEList(counts = t(counts))
  design <- matrix(1, nrow = ncol(y), ncol = 1L)
  filterByExpr(y, design = design)
}

fit_voom <- function(counts, donors, batch = FALSE) {
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
  voom <- voomWithQualityWeights(y, design, plot = FALSE)
  fit <- eBayes(lmFit(voom, design), robust = TRUE)
  table <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  answer <- data.table(
    peak_coordinate = rownames(table),
    logFC = table$logFC,
    standard_error = fit$stdev.unscaled[, coefficient] * fit$sigma,
    pvalue = table$P.Value,
    qvalue = p.adjust(table$P.Value, method = "BH")
  )
  quality_weights <- if ("sample.weights" %in% names(voom$targets)) {
    voom$targets$sample.weights
  } else {
    rep(NA_real_, nrow(donors))
  }
  attr(answer, "model_qc") <- data.table(
    donor_id = donors$donor_id,
    condition = donors$condition,
    raw_library_size = y$samples$lib.size,
    tmm_normalization_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors,
    sample_quality_weight = quality_weights,
    design_rank = qr(design)$rank,
    n_model_columns = ncol(design),
    n_tested_peaks = nrow(answer)
  )
  answer
}

fit_simple_effect <- function(counts, donors) {
  group <- factor(donors$condition, levels = c("NORMAL", "MASH"))
  design <- model.matrix(~ group)
  if (qr(design)$rank != ncol(design)) stop("Sensitivity design is not full rank")
  y <- calcNormFactors(DGEList(counts = t(counts)), method = "TMM")
  voom <- voom(y, design, plot = FALSE)
  fit <- lmFit(voom, design)
  setNames(as.numeric(fit$coefficients[, 2L]), rownames(fit$coefficients))
}

parse_peaks <- function(values) {
  match <- regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", values)
  pieces <- regmatches(values, match)
  if (any(lengths(pieces) != 4L)) stop("Malformed peak coordinate")
  data.table(
    peak_coordinate = values,
    chrom = vapply(pieces, `[`, character(1), 2L),
    start0 = as.integer(vapply(pieces, `[`, character(1), 3L)),
    end = as.integer(vapply(pieces, `[`, character(1), 4L))
  )
}

message("[DA] Importing GENCODE v49 promoters")
genes <- import(GTF, format = "gtf", feature.type = "gene")
genes <- genes[as.character(seqnames(genes)) %chin% paste0("chr", c(1:22, "X", "Y"))]
names(genes) <- as.character(mcols(genes)$gene_name)
genes <- genes[!is.na(names(genes)) & names(genes) != ""]
promoters_v49 <- promoters(genes, upstream = 2000L, downstream = 2001L)

annotate_promoters <- function(peaks) {
  parsed <- parse_peaks(peaks)
  ranges <- GRanges(
    seqnames = parsed$chrom,
    ranges = IRanges(start = parsed$start0 + 1L, end = parsed$end)
  )
  hits <- findOverlaps(ranges, promoters_v49, ignore.strand = TRUE)
  by_peak <- split(names(promoters_v49)[subjectHits(hits)], queryHits(hits))
  annotation <- rep("", length(peaks))
  for (index in names(by_peak)) {
    annotation[as.integer(index)] <- paste(sort(unique(by_peak[[index]])), collapse = ";")
  }
  annotation
}

g244_metadata <- fread(G244_META, select = c("donor_id", "batch_mm"))
joined_results <- list()
lineage_summary <- list()
model_qc_rows <- list()
filter_qc_rows <- list()

for (lineage in LINEAGES) {
  message("[DA] ", lineage)
  raw244 <- load_counts("GSE244832", lineage)
  raw281 <- load_counts("GSE281367", lineage)
  x244 <- eligible_primary(raw244)
  x281 <- eligible_primary(raw281)
  if (x244$n_normal < 4L || x244$n_mash < 4L || x281$n_normal < 4L || x281$n_mash < 4L) {
    stop("Primary lineage unexpectedly untestable: ", lineage)
  }
  keep244 <- condition_blind_filter(x244$counts)
  keep281 <- condition_blind_filter(x281$counts)
  common <- intersect(colnames(x244$counts)[keep244], colnames(x281$counts)[keep281])
  filter_qc_rows[[lineage]] <- data.table(
    release_id = RELEASE_ID,
    lineage = lineage,
    filter_method = "edgeR_filterByExpr_default_intercept_only",
    n_filter_retained_gse244832 = sum(keep244),
    n_filter_retained_gse281367 = sum(keep281),
    n_jointly_testable = length(common),
    n_normal_gse244832 = x244$n_normal,
    n_mash_gse244832 = x244$n_mash,
    n_normal_gse281367 = x281$n_normal,
    n_mash_gse281367 = x281$n_mash
  )
  if (length(common) == 0L) {
    empty <- data.table(
      peak_coordinate = character(), logFC = numeric(), standard_error = numeric(),
      pvalue = numeric(), qvalue = numeric()
    )
    fwrite(empty, file.path(OUT, paste0("cohort_primary_GSE244832_", lineage, ".tsv")), sep = "\t")
    fwrite(empty, file.path(OUT, paste0("cohort_primary_GSE281367_", lineage, ".tsv")), sep = "\t")
    lineage_summary[[lineage]] <- data.table(
      release_id = RELEASE_ID,
      lineage = lineage,
      lineage_state = "untestable",
      untestable_reason = "no_jointly_testable_peaks_after_default_filterByExpr",
      n_normal_gse244832 = x244$n_normal,
      n_mash_gse244832 = x244$n_mash,
      n_normal_gse281367 = x281$n_normal,
      n_mash_gse281367 = x281$n_mash,
      n_jointly_testable = 0L,
      n_sig_gse244832 = 0L,
      n_sig_gse281367 = 0L,
      n_supported = 0L,
      n_discordant = 0L,
      pearson_effect_correlation = NA_real_,
      spearman_effect_correlation = NA_real_,
      n_gse281367_discovery = 0L,
      n_direction_concordant = 0L,
      directional_binomial_pvalue = 1.0
    )
    next
  }
  counts244 <- x244$counts[, common, drop = FALSE]
  counts281 <- x281$counts[, common, drop = FALSE]
  fit244 <- fit_voom(counts244, x244$donors)
  fit281 <- fit_voom(counts281, x281$donors)
  qc244 <- attr(fit244, "model_qc")
  qc281 <- attr(fit281, "model_qc")
  qc244[, `:=`(release_id = RELEASE_ID, cohort = "GSE244832", lineage = lineage, model = "primary_condition_only")]
  qc281[, `:=`(release_id = RELEASE_ID, cohort = "GSE281367", lineage = lineage, model = "primary_condition_only")]
  model_qc_rows[[paste(lineage, "244", sep = "::")]] <- qc244
  model_qc_rows[[paste(lineage, "281", sep = "::")]] <- qc281
  fwrite(fit244, file.path(OUT, paste0("cohort_primary_GSE244832_", lineage, ".tsv")), sep = "\t")
  fwrite(fit281, file.path(OUT, paste0("cohort_primary_GSE281367_", lineage, ".tsv")), sep = "\t")

  joined <- merge(fit244, fit281, by = "peak_coordinate", suffixes = c("_gse244832", "_gse281367"))
  joined[, evidence_state := fifelse(
    qvalue_gse244832 < 0.05 & qvalue_gse281367 < 0.05,
    fifelse(logFC_gse244832 * logFC_gse281367 > 0, "supported", "discordant"),
    fifelse(
      xor(qvalue_gse244832 < 0.05, qvalue_gse281367 < 0.05),
      "source_dependent",
      "indeterminate"
    )
  )]
  joined[, promoter_genes := annotate_promoters(peak_coordinate)]

  # GSE244832 technical-batch sensitivity.
  batch_donors <- merge(x244$donors, g244_metadata, by = "donor_id", all.x = TRUE, sort = FALSE)
  batch_donors <- batch_donors[match(rownames(counts244), donor_id)]
  batch_fit <- fit_voom(counts244, batch_donors, batch = TRUE)
  joined[batch_fit, on = "peak_coordinate", batch_adjusted_logFC_gse244832 := i.logFC]

  # Broad GSE244832 MASL+MASH sensitivity. This never changes primary labels.
  broad_keep <- raw244$donors$contrast_eligible == TRUE & raw244$donors$condition %chin% c("NORMAL", "MASL", "MASH")
  broad_counts <- raw244$counts[broad_keep, common, drop = FALSE]
  broad_donors <- copy(raw244$donors[broad_keep])
  broad_donors[, condition := fifelse(condition == "NORMAL", "NORMAL", "MASH")]
  broad_effect <- fit_simple_effect(broad_counts, broad_donors)
  joined[, broad_masld_logFC_gse244832 := broad_effect[peak_coordinate]]

  # Leave-one-donor-out signs are evaluated only for peaks significant in either cohort.
  candidate_peaks <- joined[
    qvalue_gse244832 < 0.05 | qvalue_gse281367 < 0.05,
    peak_coordinate
  ]
  joined[, `:=`(lodo_sign_fraction_gse244832 = NA_real_, lodo_sign_fraction_gse281367 = NA_real_)]
  if (length(candidate_peaks) > 0L) {
    for (cohort in c("gse244832", "gse281367")) {
      counts <- if (cohort == "gse244832") counts244 else counts281
      donors <- if (cohort == "gse244832") x244$donors else x281$donors
      primary <- if (cohort == "gse244832") fit244 else fit281
      sign_matrix <- matrix(NA_real_, nrow = length(candidate_peaks), ncol = nrow(donors))
      rownames(sign_matrix) <- candidate_peaks
      for (drop_index in seq_len(nrow(donors))) {
        keep <- seq_len(nrow(donors)) != drop_index
        remaining <- donors[keep, .N, by = condition]
        if (any(remaining$N < 3L)) next
        effect <- fit_simple_effect(counts[keep, candidate_peaks, drop = FALSE], donors[keep])
        sign_matrix[, drop_index] <- sign(effect[candidate_peaks])
      }
      primary_sign <- sign(primary$logFC[match(candidate_peaks, primary$peak_coordinate)])
      fractions <- rowMeans(sign_matrix == primary_sign, na.rm = TRUE)
      field <- paste0("lodo_sign_fraction_", cohort)
      joined[match(candidate_peaks, peak_coordinate), (field) := fractions]
    }
  }
  joined[, sensitivity_no_reversal := (
    is.na(batch_adjusted_logFC_gse244832) |
      sign(batch_adjusted_logFC_gse244832) == sign(logFC_gse244832)
  ) & (
    is.na(lodo_sign_fraction_gse244832) | lodo_sign_fraction_gse244832 >= 0.80
  ) & (
    is.na(lodo_sign_fraction_gse281367) | lodo_sign_fraction_gse281367 >= 0.80
  )]
  joined[, lineage := lineage]
  joined_results[[lineage]] <- joined

  discovery <- joined[qvalue_gse281367 < 0.05]
  n_discovery <- nrow(discovery)
  n_same <- sum(sign(discovery$logFC_gse244832) == sign(discovery$logFC_gse281367))
  binom_p <- if (n_discovery > 0L) {
    binom.test(n_same, n_discovery, p = 0.5, alternative = "greater")$p.value
  } else 1.0
  lineage_summary[[lineage]] <- data.table(
    release_id = RELEASE_ID,
    lineage = lineage,
    lineage_state = "testable",
    untestable_reason = "",
    n_normal_gse244832 = x244$n_normal,
    n_mash_gse244832 = x244$n_mash,
    n_normal_gse281367 = x281$n_normal,
    n_mash_gse281367 = x281$n_mash,
    n_jointly_testable = nrow(joined),
    n_sig_gse244832 = sum(joined$qvalue_gse244832 < 0.05),
    n_sig_gse281367 = sum(joined$qvalue_gse281367 < 0.05),
    n_supported = sum(joined$evidence_state == "supported"),
    n_discordant = sum(joined$evidence_state == "discordant"),
    pearson_effect_correlation = cor(joined$logFC_gse244832, joined$logFC_gse281367, method = "pearson"),
    spearman_effect_correlation = cor(joined$logFC_gse244832, joined$logFC_gse281367, method = "spearman"),
    n_gse281367_discovery = n_discovery,
    n_direction_concordant = n_same,
    directional_binomial_pvalue = binom_p
  )
}

all_joined <- rbindlist(joined_results, fill = TRUE)
setcolorder(all_joined, c("lineage", "peak_coordinate", setdiff(names(all_joined), c("lineage", "peak_coordinate"))))
fwrite(all_joined, file.path(OUT, "da_peak_results.tsv.gz"), sep = "\t")
fwrite(rbindlist(filter_qc_rows), file.path(OUT, "primary_peak_filter_qc.tsv"), sep = "\t")

summary <- rbindlist(lineage_summary)
summary[, directional_binomial_qvalue := p.adjust(directional_binomial_pvalue, method = "BH")]
stability <- all_joined[, .(
  supported_sensitivity_stable = any(evidence_state == "supported") &&
    all(sensitivity_no_reversal[evidence_state == "supported"])
), by = lineage]
summary[stability, on = "lineage", supported_sensitivity_stable := i.supported_sensitivity_stable]
summary[is.na(supported_sensitivity_stable), supported_sensitivity_stable := FALSE]
summary[, dynamic_main_eligible := n_supported > 0L & directional_binomial_qvalue < 0.05 &
  supported_sensitivity_stable]
fwrite(summary, file.path(OUT, "da_lineage_summary.tsv"), sep = "\t")
model_qc <- rbindlist(model_qc_rows, fill = TRUE)
setcolorder(model_qc, c(
  "release_id", "cohort", "lineage", "model", "donor_id", "condition",
  "raw_library_size", "tmm_normalization_factor", "effective_library_size",
  "sample_quality_weight", "design_rank", "n_model_columns", "n_tested_peaks"
))
fwrite(model_qc, file.path(OUT, "da_model_qc.tsv"), sep = "\t")

# Historical GSE244832 peak coordinates are a sensitivity only. They never alter
# the primary jointly testable family or its evidence states.
legacy_root <- file.path(CANDIDATE, "counts_legacy_gse244832")
if (dir.exists(legacy_root)) {
  legacy_rows <- list()
  legacy_peak_rows <- list()
  for (lineage in LINEAGES) {
    legacy244 <- eligible_primary(load_counts("GSE244832", lineage, root_override = legacy_root))
    legacy281 <- eligible_primary(load_counts("GSE281367", lineage, root_override = legacy_root))
    keep244 <- condition_blind_filter(legacy244$counts)
    keep281 <- condition_blind_filter(legacy281$counts)
    common <- intersect(colnames(legacy244$counts)[keep244], colnames(legacy281$counts)[keep281])
    if (length(common) == 0L) {
      legacy_rows[[lineage]] <- data.table(
        release_id = RELEASE_ID,
        lineage = lineage,
        state = "untestable",
        reason = "no_jointly_testable_peaks_after_default_filterByExpr",
        n_jointly_testable = 0L,
        pearson_effect_correlation = NA_real_,
        spearman_effect_correlation = NA_real_,
        n_same_direction = 0L,
        n_supported = 0L
      )
      next
    }
    a <- fit_voom(legacy244$counts[, common, drop = FALSE], legacy244$donors)
    b <- fit_voom(legacy281$counts[, common, drop = FALSE], legacy281$donors)
    joined <- merge(a, b, by = "peak_coordinate", suffixes = c("_gse244832", "_gse281367"))
    joined[, `:=`(
      release_id = RELEASE_ID,
      lineage = lineage,
      evidence_state = fcase(
        qvalue_gse244832 < 0.05 & qvalue_gse281367 < 0.05 &
          logFC_gse244832 * logFC_gse281367 > 0, "supported",
        qvalue_gse244832 < 0.05 & qvalue_gse281367 < 0.05 &
          logFC_gse244832 * logFC_gse281367 < 0, "discordant",
        xor(qvalue_gse244832 < 0.05, qvalue_gse281367 < 0.05), "source_dependent",
        default = "indeterminate"
      )
    )]
    legacy_peak_rows[[lineage]] <- joined
    legacy_rows[[lineage]] <- data.table(
      release_id = RELEASE_ID,
      lineage = lineage,
      state = "testable",
      reason = "",
      n_jointly_testable = nrow(joined),
      pearson_effect_correlation = cor(joined$logFC_gse244832, joined$logFC_gse281367),
      spearman_effect_correlation = cor(joined$logFC_gse244832, joined$logFC_gse281367, method = "spearman"),
      n_same_direction = sum(sign(joined$logFC_gse244832) == sign(joined$logFC_gse281367)),
      n_supported = sum(
        joined$qvalue_gse244832 < 0.05 & joined$qvalue_gse281367 < 0.05 &
          joined$logFC_gse244832 * joined$logFC_gse281367 > 0
      )
    )
  }
  fwrite(rbindlist(legacy_rows), file.path(OUT, "legacy_gse244832_peak_space_sensitivity.tsv"), sep = "\t")
  fwrite(rbindlist(legacy_peak_rows, fill = TRUE), file.path(OUT, "legacy_gse244832_peak_space_peak_results.tsv.gz"), sep = "\t")
} else {
  fwrite(data.table(
    release_id = RELEASE_ID,
    sensitivity = "legacy_gse244832_peak_space",
    state = "untestable",
    reason = "legacy_exact_counts_unavailable"
  ), file.path(OUT, "legacy_gse244832_peak_space_sensitivity.tsv"), sep = "\t")
}

# Macrophage 50-cell sensitivity is tested only if both cohorts retain four
# donors per condition; otherwise its explicit underpowered state is retained.
mac244 <- load_counts("GSE244832", "macrophage")
mac281 <- load_counts("GSE281367", "macrophage")
strict_counts <- function(value, target_condition) {
  answer <- value$donors[n_cells >= 50 & condition == target_condition, .N]
  if (length(answer) == 0L || is.na(answer)) 0L else answer
}
strict_244_normal <- strict_counts(mac244, "NORMAL")
strict_244_mash <- strict_counts(mac244, "MASH")
strict_281_normal <- strict_counts(mac281, "NORMAL")
strict_281_mash <- strict_counts(mac281, "MASH")
strict_testable <- min(
  strict_244_normal, strict_244_mash, strict_281_normal, strict_281_mash
) >= 4L
fwrite(
  data.table(
    release_id = RELEASE_ID,
    lineage = "macrophage",
    min_cells = 50L,
    n_mash_gse244832 = strict_244_mash,
    n_normal_gse244832 = strict_244_normal,
    n_mash_gse281367 = strict_281_mash,
    n_normal_gse281367 = strict_281_normal,
    testable = strict_testable,
    state = ifelse(strict_testable, "testable", "untestable"),
    reason = ifelse(strict_testable, "", "fewer_than_four_donors_per_condition_in_at_least_one_cohort")
  ),
  file.path(OUT, "macrophage_50cell_sensitivity.tsv"),
  sep = "\t"
)
if (strict_testable) {
  keep244 <- mac244$donors$n_cells >= 50 & mac244$donors$condition %chin% c("NORMAL", "MASH")
  keep281 <- mac281$donors$n_cells >= 50 & mac281$donors$condition %chin% c("NORMAL", "MASH")
  counts244 <- mac244$counts[keep244, , drop = FALSE]
  counts281 <- mac281$counts[keep281, , drop = FALSE]
  common <- intersect(
    colnames(counts244)[condition_blind_filter(counts244)],
    colnames(counts281)[condition_blind_filter(counts281)]
  )
  strict_a <- fit_voom(counts244[, common, drop = FALSE], mac244$donors[keep244])
  strict_b <- fit_voom(counts281[, common, drop = FALSE], mac281$donors[keep281])
  strict_joined <- merge(strict_a, strict_b, by = "peak_coordinate", suffixes = c("_gse244832", "_gse281367"))
  fwrite(strict_joined, file.path(OUT, "macrophage_50cell_peak_results.tsv.gz"), sep = "\t")
}

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("[DA] Wrote ", OUT)
