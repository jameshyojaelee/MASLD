#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

LNCRNA_COHORTS <- c(
  "GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"
)
LNCRNA_TREAT_LFC <- 0.25   # retained: interval-null comparator arm only
LNCRNA_FDR <- 0.05
# CANONICAL 2026-08-12: padj < 0.05 AND |log2FC| > 0.50, one floor for every
# biotype. A biotype-specific floor was considered and rejected: lncRNAs are
# lower-expressed, so their logFC estimates are noisier, and a hard floor on the
# point estimate already admits them ABOVE their share of the tested universe
# (29.8% of positives vs 26.7% of genes at 0.50, rising to 32.1% at 0.20).
# Loosening the lncRNA floor would amplify that, and a class-specific filter
# would make the biotype comparison circular by construction.
LNCRNA_LFC <- 0.50

# Canonical positive call. Accepts either adjusted-p column name: the pooled
# primary table uses `padj`, fit_disease_contrast() emits `FDR`.
lncrna_canonical_positive <- function(dt, fdr = LNCRNA_FDR, lfc = LNCRNA_LFC) {
  p_col <- intersect(c("padj", "FDR", "adj.P.Val"), names(dt))[1L]
  if (is.na(p_col)) stop("lncrna_canonical_positive: no adjusted p-value column", call. = FALSE)
  !is.na(dt[[p_col]]) & !is.na(dt$logFC) & dt[[p_col]] < fdr & abs(dt$logFC) > lfc
}

fail <- function(...) stop(..., call. = FALSE)

is_symlink <- function(path) {
  target <- Sys.readlink(path)
  !is.na(target) && nzchar(target)
}

assert_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
}

read_source_gate <- function(path) {
  if (!file.exists(path) || is_symlink(path)) {
    fail("Missing or symlinked source gate: ", path)
  }
  gate <- fread(path, colClasses = "character")
  assert_true(
    identical(names(gate), c("gate", "status", "detail")),
    "Unexpected source-gate schema"
  )
  assert_true(!anyDuplicated(gate$gate), "Duplicate source-gate rows")
  selector <- gate[["gate"]] == "lncrna_source_gate"
  row <- gate[selector]
  assert_true(nrow(row) == 1L, "Missing unique lncrna_source_gate row")
  assert_true(
    identical(row$status, "pass"),
    paste0("lncrna_source_gate is not pass: ", row$status)
  )
  invisible(gate)
}

make_disease_design <- function(samples, include_dataset) {
  required <- c("sample_id", "dataset", "inferred_sex", "group_binary")
  assert_true(all(required %in% names(samples)), "Model metadata fields missing")
  assert_true(!anyDuplicated(samples$sample_id), "Duplicate model sample IDs")
  samples <- copy(samples)
  samples[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]
  samples[, inferred_sex := droplevels(factor(inferred_sex))]
  assert_true(
    !anyNA(samples$group_binary) && nlevels(samples$group_binary) == 2L,
    "Disease contrast is not estimable"
  )
  terms <- character()
  if (include_dataset) {
    samples[, dataset := droplevels(factor(dataset))]
    assert_true(nlevels(samples$dataset) >= 2L, "Pooled design lacks multiple cohorts")
    terms <- c(terms, "dataset")
  }
  if (nlevels(samples$inferred_sex) >= 2L) terms <- c(terms, "inferred_sex")
  terms <- c(terms, "group_binary")
  design <- model.matrix(
    as.formula(paste("~", paste(terms, collapse = " + "))),
    data = samples
  )
  rownames(design) <- samples$sample_id
  assert_true(
    "group_binaryDisease" %in% colnames(design),
    "Disease coefficient missing"
  )
  assert_true(qr(design)$rank == ncol(design), "Model design is rank deficient")
  design
}

fit_disease_contrast <- function(
    dge,
    samples,
    include_dataset,
    quality_weights = TRUE,
    treat_lfc = LNCRNA_TREAT_LFC) {
  assert_true(identical(colnames(dge), samples$sample_id), "DGE/metadata order drift")
  design <- make_disease_design(samples, include_dataset)
  if (quality_weights) {
    transformed <- voomWithQualityWeights(dge, design, plot = FALSE)
  } else {
    transformed <- voom(dge, design, plot = FALSE)
  }
  fit0 <- lmFit(transformed, design)
  coefficient <- match("group_binaryDisease", colnames(design))
  standard_fit <- eBayes(fit0)
  standard_error <- standard_fit$stdev.unscaled[, coefficient] *
    sqrt(standard_fit$s2.post)
  critical <- qt(0.975, df = standard_fit$df.total)
  names(standard_error) <- rownames(standard_fit$coefficients)
  names(critical) <- rownames(standard_fit$coefficients)
  standard <- topTable(
    standard_fit,
    coef = coefficient,
    number = Inf,
    sort.by = "none"
  )
  treat_result <- topTreat(
    treat(fit0, lfc = treat_lfc),
    coef = coefficient,
    number = Inf,
    sort.by = "none"
  )
  genes <- rownames(dge)
  result <- data.table(
    gene_id_versioned = genes,
    logFC = standard[genes, "logFC"],
    SE = standard_error[genes],
    CI_low = standard[genes, "logFC"] - critical[genes] * standard_error[genes],
    CI_high = standard[genes, "logFC"] + critical[genes] * standard_error[genes],
    P.Value = standard[genes, "P.Value"],
    FDR = standard[genes, "adj.P.Val"],
    treat_lfc = treat_lfc,
    treat_p = treat_result[genes, "P.Value"],
    treat_fdr = treat_result[genes, "adj.P.Val"],
    AveExpr = standard[genes, "AveExpr"]
  )
  numeric_fields <- setdiff(names(result), "gene_id_versioned")
  assert_true(
    nrow(result) == nrow(dge) && !anyDuplicated(result$gene_id_versioned),
    "Model output gene-universe drift"
  )
  assert_true(
    all(vapply(numeric_fields, function(field) {
      all(is.finite(result[[field]]))
    }, logical(1))),
    "Model output contains non-finite values"
  )
  result
}

direction_code <- function(value) {
  fifelse(value > 0, 1L, fifelse(value < 0, -1L, 0L))
}

derive_high_confidence <- function(primary, cohort_long, loco_long, equal_weight) {
  required_primary <- c(
    "gene_id_versioned", "logFC", "padj", "treat_fdr", "is_canonical_chromosome",
    "mapping_status", "main_text_eligible"
  )
  assert_true(all(required_primary %in% names(primary)), "Primary high-confidence fields missing")
  assert_true(!anyDuplicated(primary$gene_id_versioned), "Duplicate primary genes")
  assert_true(
    uniqueN(cohort_long$cohort) == 5L &&
      all(cohort_long[, .N, by = gene_id_versioned]$N == 5L),
    "Per-cohort family is incomplete"
  )
  assert_true(
    uniqueN(loco_long$excluded_cohort) == 5L &&
      all(loco_long[, .N, by = gene_id_versioned]$N == 5L),
    "Leave-one-cohort-out family is incomplete"
  )
  assert_true(!anyDuplicated(equal_weight$gene_id_versioned), "Duplicate equal-weight genes")

  base <- copy(primary)
  base[, pooled_direction := direction_code(logFC)]
  cohort <- merge(
    cohort_long,
    base[, .(gene_id_versioned, pooled_direction)],
    by = "gene_id_versioned",
    all.x = TRUE
  )
  cohort[, `:=`(
    direction_agrees = direction_code(logFC) == pooled_direction,
    materially_opposite = fifelse(
      pooled_direction > 0L,
      CI_high < 0,
      fifelse(pooled_direction < 0L, CI_low > 0, TRUE)
    )
  )]
  cohort_summary <- cohort[, .(
    cohort_direction_agreement_n = sum(direction_agrees),
    cohort_materially_opposite_n = sum(materially_opposite),
    cohort_all_finite = all(is.finite(logFC) & is.finite(CI_low) & is.finite(CI_high))
  ), by = gene_id_versioned]

  loco <- merge(
    loco_long,
    base[, .(gene_id_versioned, pooled_direction)],
    by = "gene_id_versioned",
    all.x = TRUE
  )
  loco[, direction_agrees := direction_code(logFC) == pooled_direction]
  loco_summary <- loco[, .(
    loco_direction_agreement_n = sum(direction_agrees),
    loco_canonical_positive_n = sum(FDR < LNCRNA_FDR & abs(logFC) > LNCRNA_LFC),
    loco_treat_fdr_lt_0.05_n = sum(treat_fdr < LNCRNA_FDR)   # comparator arm
  ), by = gene_id_versioned]

  equal <- merge(
    equal_weight[, .(
      gene_id_versioned,
      equal_weight_logFC = logFC,
      equal_weight_FDR = FDR,
      equal_weight_treat_p = treat_p,
      equal_weight_treat_fdr = treat_fdr
    )],
    base[, .(gene_id_versioned, pooled_direction)],
    by = "gene_id_versioned"
  )
  equal[, equal_weight_direction_agrees :=
          direction_code(equal_weight_logFC) == pooled_direction]

  result <- Reduce(
    function(left, right) merge(left, right, by = "gene_id_versioned", all.x = TRUE),
    list(
      base,
      cohort_summary,
      loco_summary,
      equal[, setdiff(names(equal), "pooled_direction"), with = FALSE]
    )
  )
  result[, primary_canonical_positive := lncrna_canonical_positive(result)]
  result[, primary_treat_positive := treat_fdr < LNCRNA_FDR]   # comparator arm
  result[, genomic_mapping_eligible :=
           is_canonical_chromosome == "true" &
           mapping_status == "mapping_unambiguous" &
           main_text_eligible == "true"]
  result[, high_confidence :=
           primary_canonical_positive &
           genomic_mapping_eligible &
           cohort_direction_agreement_n >= 4L &
           cohort_materially_opposite_n == 0L &
           cohort_all_finite &
           loco_direction_agreement_n == 5L &
           loco_canonical_positive_n >= 4L &
           equal_weight_direction_agrees &
           equal_weight_FDR < LNCRNA_FDR &
           abs(equal_weight_logFC) > LNCRNA_LFC]
  result[]
}
