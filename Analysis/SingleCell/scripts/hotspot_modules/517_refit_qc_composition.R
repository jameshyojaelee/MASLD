#!/usr/bin/env Rscript

# Candidate-only Figure 2B correction. Refit the existing MuSiC proportion
# analysis after enforcing the canonical pass_technical gate. This appends a
# Plan 20 artifact; it never writes the live celltype_attribution results.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) {
    stop("Missing required argument: ", flag)
  }
  args[[position + 1L]]
}

project_root <- normalizePath(arg_value("--project-root"), mustWork = TRUE)
output_root <- normalizePath(arg_value("--output-root"), mustWork = TRUE)
producer_script <- normalizePath(arg_value("--producer-script"), mustWork = TRUE)
candidate_id <- "program-context-v2-candidate-2026-08-07"

proportions_path <- file.path(
  project_root,
  "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"
)
qc_path <- file.path(
  project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"
)
legacy_path <- file.path(
  project_root,
  "RNA-seq/results/celltype_attribution/composition_shifts.csv"
)
input_paths <- c(proportions_path, qc_path, legacy_path, producer_script)
if (any(!file.exists(input_paths))) {
  stop("A required composition input is missing")
}

output_paths <- file.path(
  output_root,
  c(
    "composition_sample_qc_v2.tsv",
    "composition_sample_qc_testability.tsv",
    "composition_sample_qc_sensitivity.tsv",
    "composition_sample_qc_input_manifest.tsv",
    "composition_sample_qc_audit.tsv"
  )
)
if (any(file.exists(output_paths))) {
  stop("Refusing to overwrite an existing Plan 20 composition artifact")
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path)
  }
  sub(" .*", "", output[[1L]])
}

atomic_fwrite <- function(table, path) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  fwrite(table, temporary, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    stop("Atomic rename failed for ", path)
  }
}

props <- fread(proportions_path, check.names = FALSE)
qc <- fread(qc_path, check.names = FALSE)
key <- c("sample_id", "dataset")
if (anyDuplicated(props, by = key) || anyDuplicated(qc, by = key)) {
  stop("Composition or QC input has duplicate sample_id/dataset keys")
}
required_props <- c(key, "group_binary", "inferred_sex")
required_qc <- c(key, "pass_technical")
if (!all(required_props %in% names(props)) || !all(required_qc %in% names(qc))) {
  stop("Composition or QC input lacks required columns")
}

qc[, pass_technical := as.logical(pass_technical)]
eligible <- qc[pass_technical == TRUE]
if (nrow(eligible) != 1260L || uniqueN(eligible$dataset) != 9L) {
  stop("Canonical pass_technical census drift: expected 1,260 samples/9 cohorts")
}

joined_all <- merge(
  props,
  qc[, c(key, "pass_technical"), with = FALSE],
  by = key,
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(joined_all$pass_technical)) {
  stop("A deconvolved sample is absent from the canonical QC report")
}
qc_failed_excluded <- joined_all[pass_technical == FALSE]
analyzed <- joined_all[pass_technical == TRUE]
missing_deconvolution <- fsetdiff(
  eligible[, ..key],
  props[, ..key],
  all = TRUE
)
if (
  nrow(analyzed) != 1221L ||
  nrow(missing_deconvolution) != 39L ||
  nrow(qc_failed_excluded) != 21L
) {
  stop(
    "Composition/QC join drift: expected 1,221 analyzed, 39 assay-missing, ",
    "and 21 technical-QC failures excluded"
  )
}
if (!setequal(unique(analyzed$group_binary), c("Control", "Disease"))) {
  stop("QC-filtered composition data lack the Control/Disease contrast")
}

metadata_columns <- c(
  "sample_id", "dataset", "group_binary", "condition", "fibrosis_stage",
  "inferred_sex", "diagnosis_harmonized", "f2_group", "pass_technical"
)
celltypes <- setdiff(names(analyzed), metadata_columns)
if (length(celltypes) != 22L) {
  stop("Expected exactly 22 deposited cell-type columns")
}

transform_proportion <- function(value) {
  asin(sqrt(pmax(pmin(value, 1 - 1e-6), 1e-6)))
}

fit_celltype <- function(celltype) {
  values <- analyzed[[celltype]]
  keep <- !is.na(values)
  if (sum(keep) < 30L) {
    return(NULL)
  }
  data <- analyzed[keep]
  outcome <- transform_proportion(data[[celltype]])
  group <- factor(data$group_binary, levels = c("Control", "Disease"))
  dataset <- factor(data$dataset)
  inferred_sex <- factor(data$inferred_sex)
  if (anyNA(inferred_sex)) {
    stop("Missing inferred_sex in a testable composition model: ", celltype)
  }
  design <- model.matrix(~ group + dataset + inferred_sex)
  if (nrow(design) != length(outcome) || !"groupDisease" %in% colnames(design)) {
    stop("Composition design failure: ", celltype)
  }
  fit <- eBayes(lmFit(matrix(outcome, nrow = 1L), design))
  result <- topTable(fit, coef = "groupDisease", number = 1L, sort.by = "none")
  data.table(
    celltype = celltype,
    testability = "testable",
    n_total_qc = nrow(eligible),
    n_analyzed = length(outcome),
    n_control = sum(group == "Control"),
    n_disease = sum(group == "Disease"),
    effect_arcsin_sqrt = result$logFC,
    t = result$t,
    pvalue = result$P.Value
  )
}

fit_rows <- Filter(Negate(is.null), lapply(celltypes, fit_celltype))
results <- rbindlist(fit_rows, use.names = TRUE, fill = TRUE)
if (nrow(results) != 16L || anyDuplicated(results$celltype)) {
  stop("Expected a complete 16-cell-type testable family")
}
results[, padj := p.adjust(pvalue, method = "BH")]
results[, `:=`(
  model_formula = "asin(sqrt(proportion)) ~ group_binary + dataset + inferred_sex",
  biological_unit = "QC-passing human liver sample with MuSiC deconvolution",
  multiplicity_family = "BH across all 16 testable deposited cell-type columns",
  source_assay = "bulk RNA-seq MuSiC deconvolution",
  allowed_wording = "QC-filtered sample-level cross-sectional composition association",
  prohibited_wording = "donor-level, longitudinal change, lineage transition, or causal composition effect"
)]
setorder(results, celltype)

testability <- data.table(celltype = celltypes)
testability[, n_nonmissing_qc := vapply(
  celltype,
  function(celltype) sum(!is.na(analyzed[[celltype]])),
  integer(1L)
)]
testability[, testability := fifelse(
  n_nonmissing_qc >= 30L,
  "testable",
  "untestable_structurally_unavailable"
)]
testability[, reason := fifelse(
  testability == "testable",
  "included in the complete BH family",
  "deposited merged column is structurally all-NA"
)]
setorder(testability, celltype)
if (
  testability[testability == "testable", .N] != 16L ||
  testability[testability != "testable", .N] != 6L
) {
  stop("Cell-type testability census drift")
}

legacy <- fread(legacy_path)[contrast == "Control_vs_Disease"]
sensitivity <- merge(
  results[, .(
    celltype,
    refit_effect = effect_arcsin_sqrt,
    refit_pvalue = pvalue,
    refit_padj = padj
  )],
  legacy[, .(
    celltype,
    legacy_effect = logit_diff,
    legacy_pvalue = pvalue,
    legacy_padj = padj
  )],
  by = "celltype",
  all = TRUE
)
if (nrow(sensitivity) != 16L || anyNA(sensitivity$refit_effect) || anyNA(sensitivity$legacy_effect)) {
  stop("Legacy/refit sensitivity join is incomplete")
}
sensitivity[, `:=`(
  direction_agreement = sign(refit_effect) == sign(legacy_effect),
  fdr_call_agreement = (refit_padj < 0.05) == (legacy_padj < 0.05)
)]
setorder(sensitivity, celltype)

input_manifest <- data.table(
  artifact_id = c(
    "composition_proportions", "canonical_qc_report",
    "legacy_composition_sensitivity_only", "producer_script"
  ),
  source_path = substring(input_paths, nchar(project_root) + 2L),
  sha256 = vapply(input_paths, sha256_file, character(1L)),
  bytes = file.info(input_paths)$size,
  role = c(
    "analysis_input", "analysis_input", "sensitivity_only", "producer"
  )
)

audit <- data.table(
  audit_id = c(
    "pass_technical_census", "qc_join", "deconvolution_missingness",
    "excluded_qc_failures", "declared_celltypes", "testable_family",
    "untestable_columns", "bh_family", "direction_sensitivity",
    "fdr_sensitivity", "biological_unit", "canonical_write"
  ),
  status = "pass",
  value = c(
    "1260", "1221", "39", "21", "22", "16", "6", "16",
    as.character(sum(sensitivity$direction_agreement)),
    as.character(sum(sensitivity$fdr_call_agreement)),
    "human liver sample", "false"
  ),
  detail = c(
    "Canonical pass_technical samples across 9 retained cohorts",
    "QC-passing samples with deposited MuSiC proportions",
    "QC-passing samples without deposited deconvolution",
    "Deconvolved samples excluded for failing pass_technical",
    "All deposited cell-type columns enumerated before outcome testing",
    "Complete inferential family with non-degenerate deposited values",
    "Structurally all-NA deposited columns retained as untestable",
    "BH applied once across all 16 testable cell types",
    "Refit/legacy effect-sign agreement out of 16",
    "Refit/legacy FDR-call agreement out of 16",
    "No authoritative donor identifier was used; donor-level wording prohibited",
    "Candidate-only Plan 20 append; no canonical result overwritten"
  )
)

atomic_fwrite(results, output_paths[[1L]])
atomic_fwrite(testability, output_paths[[2L]])
atomic_fwrite(sensitivity, output_paths[[3L]])
atomic_fwrite(input_manifest, output_paths[[4L]])
atomic_fwrite(audit, output_paths[[5L]])

message("Wrote candidate-only QC-filtered composition artifacts to ", output_root)
