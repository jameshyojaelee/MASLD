#!/usr/bin/env Rscript

# Independent candidate-only validation for 517_refit_qc_composition.R.
# Re-derives sample/testability/BH counts and checks each limma coefficient
# against a base-R OLS fit with the identical design.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) stop("Missing argument: ", flag)
  args[[position + 1L]]
}

project_root <- normalizePath(arg_value("--project-root"), mustWork = TRUE)
output_root <- normalizePath(arg_value("--output-root"), mustWork = TRUE)
producer_script <- normalizePath(arg_value("--producer-script"), mustWork = TRUE)
candidate_id <- "program-context-v2-candidate-2026-08-07"

result_path <- file.path(output_root, "composition_sample_qc_v2.tsv")
testability_path <- file.path(output_root, "composition_sample_qc_testability.tsv")
audit_path <- file.path(output_root, "composition_sample_qc_audit.tsv")
validation_path <- file.path(output_root, "composition_sample_qc_validation.tsv")
seal_path <- file.path(output_root, "COMPOSITION_QC_READY")
if (any(!file.exists(c(result_path, testability_path, audit_path)))) {
  stop("Required composition artifacts are missing")
}
if (file.exists(validation_path) || file.exists(seal_path)) {
  stop("Refusing to overwrite composition validation or ready seal")
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed: ", path)
  sub(" .*", "", output[[1L]])
}
atomic_fwrite <- function(table, path) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  fwrite(table, temporary, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Atomic rename failed: ", path)
}

props <- fread(file.path(
  project_root,
  "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"
), check.names = FALSE)
qc <- fread(file.path(
  project_root,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"
), check.names = FALSE)
results <- fread(result_path)
testability <- fread(testability_path)
audit <- fread(audit_path)
key <- c("sample_id", "dataset")
qc[, pass_technical := as.logical(pass_technical)]
eligible <- qc[pass_technical == TRUE]
analyzed <- merge(
  props,
  eligible[, ..key],
  by = key,
  all = FALSE,
  sort = FALSE
)
missing <- fsetdiff(eligible[, ..key], props[, ..key], all = TRUE)
if (
  nrow(eligible) != 1260L || nrow(analyzed) != 1221L || nrow(missing) != 39L ||
  nrow(results) != 16L || nrow(testability) != 22L ||
  testability[testability == "testable", .N] != 16L
) {
  stop("Independent composition census validation failed")
}

transform_proportion <- function(value) {
  asin(sqrt(pmax(pmin(value, 1 - 1e-6), 1e-6)))
}
ols_effects <- vapply(results$celltype, function(celltype) {
  values <- analyzed[[celltype]]
  keep <- !is.na(values)
  data <- analyzed[keep]
  data[, outcome := transform_proportion(get(celltype))]
  data[, group := factor(group_binary, levels = c("Control", "Disease"))]
  fit <- lm(outcome ~ group + factor(dataset) + factor(inferred_sex), data = data)
  unname(coef(fit)[["groupDisease"]])
}, numeric(1L))
max_effect_delta <- max(abs(ols_effects - results$effect_arcsin_sqrt))
if (!is.finite(max_effect_delta) || max_effect_delta > 1e-10) {
  stop("Independent OLS coefficients do not reproduce limma effects")
}

bh_rederived <- p.adjust(results$pvalue, method = "BH")
max_bh_delta <- max(abs(bh_rederived - results$padj))
if (!is.finite(max_bh_delta) || max_bh_delta > 1e-12) {
  stop("Independent BH rederivation failed")
}
required_audits <- c(
  pass_technical_census = "1260", qc_join = "1221",
  deconvolution_missingness = "39", excluded_qc_failures = "21",
  declared_celltypes = "22", testable_family = "16",
  untestable_columns = "6", bh_family = "16",
  biological_unit = "human liver sample", canonical_write = "false"
)
observed_audits <- setNames(as.character(audit$value), audit$audit_id)
if (
  !all(names(required_audits) %in% names(observed_audits)) ||
  !identical(unname(observed_audits[names(required_audits)]), unname(required_audits))
) {
  stop("Composition audit facts are incomplete or stale")
}
if (any(grepl("donor-level", results$allowed_wording, fixed = TRUE))) {
  stop("Candidate composition table contains donor-level wording")
}

validation <- data.table(
  check_id = c(
    "qc_census", "assay_join", "declared_universe", "ols_effects",
    "bh_rederivation", "sample_level_wording", "candidate_only"
  ),
  status = "pass",
  observed = c(
    "1260", "1221_with_39_missing", "22_declared_16_testable_6_untestable",
    format(max_effect_delta, scientific = TRUE),
    format(max_bh_delta, scientific = TRUE), "human_liver_sample", "true"
  ),
  criterion = c(
    "1260 pass_technical across 9 cohorts",
    "1221 QC-pass deconvolved; 39 QC-pass assay-missing",
    "all deposited cell-type columns explicitly classified",
    "base-R OLS group coefficients agree within 1e-10",
    "BH across complete 16-test family agrees within 1e-12",
    "donor-level wording absent",
    "no live canonical output written"
  )
)
atomic_fwrite(validation, validation_path)

seal <- data.table(
  candidate_id = candidate_id,
  status = "QC_FILTERED_SAMPLE_LEVEL_COMPOSITION_READY",
  n_pass_technical = 1260L,
  n_deconvolved_qc_pass = 1221L,
  n_assay_missing = 39L,
  n_qc_failures_excluded = 21L,
  n_declared_celltypes = 22L,
  n_testable_celltypes = 16L,
  n_untestable_celltypes = 6L,
  model = "asin(sqrt(proportion)) ~ group_binary + dataset + inferred_sex",
  biological_unit = "human liver sample",
  multiplicity_family = "BH across all 16 testable deposited cell-type columns",
  results_sha256 = sha256_file(result_path),
  testability_sha256 = sha256_file(testability_path),
  audit_sha256 = sha256_file(audit_path),
  validation_sha256 = sha256_file(validation_path),
  validator_sha256 = sha256_file(producer_script),
  canonical_promotion_authorized = "false"
)
atomic_fwrite(seal, seal_path)
message("QC-filtered sample-level composition validation passed")
