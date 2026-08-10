#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: validate_arm.R RUN_ROOT ARM", call. = FALSE)
run_root <- normalizePath(args[[1]], mustWork = TRUE)
arm <- args[[2]]
arms <- c("R0", "F_locked", "F_legacy", "F_five")
if (!arm %in% arms) stop("Unknown BG-001 arm: ", arm, call. = FALSE)
if (!file.exists(file.path(run_root, ".bg001_candidate_root"))) stop("Missing candidate sentinel", call. = FALSE)

arm_root <- file.path(run_root, "arms", arm)
contract <- jsonlite::read_json(file.path(run_root, "contract/run_contract.json"), simplifyVector = FALSE)
mode <- contract$arm_modes[[arm]]
project_root <- file.path(run_root, "source_snapshot")
cfg <- yaml::read_yaml(file.path(project_root, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(x) isTRUE(x$de$include_in_mega), cfg))
sha256 <- function(path) {
  line <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(line, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path, call. = FALSE)
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

dge <- readRDS(file.path(arm_root, "results/integration/merged_dge.rds"))
raw <- readRDS(file.path(arm_root, "results/integration/merged_counts_raw.rds"))
meta <- as.data.table(readRDS(file.path(arm_root, "results/integration/meta_matched.rds")))
deg <- fread(file.path(arm_root, "results/integration/deg_results.csv"))
weights <- fread(file.path(arm_root, "results/integration/lvqw_sample_weights.tsv"))
design_observed <- fread(file.path(arm_root, "results/integration/model_design.tsv"))
model_manifest <- fread(file.path(arm_root, "results/integration/model_input_manifest.tsv"))
preprocessing <- fread(file.path(arm_root, "provenance/preprocessing_manifest.tsv"))
filter_stats <- fread(file.path(arm_root, "provenance/gene_filter_statistics.tsv.gz"))
effective_sources <- fread(file.path(arm_root, "provenance/effective_count_sources.tsv"))

if (anyDuplicated(rownames(dge)) || anyDuplicated(colnames(dge)) ||
    anyDuplicated(deg$gene) || anyDuplicated(meta$sample_id) ||
    anyDuplicated(weights$sample_id) || anyDuplicated(filter_stats$gene)) {
  stop("Duplicate gene/sample identifier in arm ", arm, call. = FALSE)
}
if (!identical(deg$gene, rownames(dge))) {
  stop("DEG gene IDs/order differ from the arm DGE", call. = FALSE)
}
if (!identical(rownames(dge$samples), colnames(dge))) {
  stop("DGE sample annotation order differs from count columns", call. = FALSE)
}
if (any(!is.finite(dge$samples$lib.size)) || any(dge$samples$lib.size <= 0) ||
    any(!is.finite(dge$samples$norm.factors)) || any(dge$samples$norm.factors <= 0)) {
  stop("DGE library sizes/norm factors are not finite positive", call. = FALSE)
}

# Validate the RNA-seq count domain in bounded row chunks to avoid allocating a
# second full logical matrix on the compute node.
for (first in seq.int(1L, nrow(dge), by = 256L)) {
  last <- min(nrow(dge), first + 255L)
  block <- dge$counts[first:last, , drop = FALSE]
  if (any(!is.finite(block)) || any(block < 0) || any(block != floor(block))) {
    stop("DGE contains a nonfinite, negative, or noninteger count", call. = FALSE)
  }
}

required_numeric <- c("logFC", "SE", "t", "P.Value", "padj", "treat_lfc", "treat_p", "treat_fdr", "AveExpr")
missing_numeric <- setdiff(required_numeric, names(deg))
if (length(missing_numeric) ||
    any(vapply(required_numeric, function(field) any(!is.finite(deg[[field]])), logical(1))) ||
    any(deg$SE <= 0) || any(deg$treat_lfc != 0.25)) {
  stop("DEG coefficient domain or exact TREAT floor is invalid", call. = FALSE)
}
for (field in c("P.Value", "padj", "treat_p", "treat_fdr")) {
  if (any(deg[[field]] < 0 | deg[[field]] > 1)) {
    stop("Probability/FDR field is outside [0,1]: ", field, call. = FALSE)
  }
}

if (nrow(preprocessing) != 1L ||
    !identical(preprocessing$counts_mode[[1L]], mode$counts) ||
    !identical(preprocessing$qc_mode[[1L]], mode$qc) ||
    !identical(preprocessing$filter_scope[[1L]], mode$filter_scope) ||
    !identical(preprocessing$gene_mode[[1L]], mode$gene_mode)) {
  stop("Preprocessing manifest differs from the frozen arm mode", call. = FALSE)
}
required_filter_fields <- c(
  "gene", "total_count", "n_samples_nonzero", "n_cohorts_nonzero", "mean_logCPM", "kept"
)
if (!all(required_filter_fields %in% names(filter_stats))) {
  stop("Filter-stat table lacks required fields", call. = FALSE)
}
if (!identical(filter_stats$gene, rownames(raw)) ||
    anyNA(filter_stats$kept) ||
    !identical(filter_stats[kept == TRUE, gene], rownames(dge))) {
  stop("Filter-stat table does not exactly cover raw genes/kept DGE genes", call. = FALSE)
}
filter_numeric <- c("total_count", "n_samples_nonzero", "n_cohorts_nonzero", "mean_logCPM")
if (any(vapply(filter_numeric, function(field) any(!is.finite(filter_stats[[field]])), logical(1))) ||
    any(filter_stats$total_count < 0) || any(filter_stats$n_samples_nonzero < 0) ||
    any(filter_stats$n_cohorts_nonzero < 0)) {
  stop("Filter-stat table has an invalid numerical domain", call. = FALSE)
}

model_samples <- colnames(dge)[dge$samples$dataset %in% mega]
model_meta <- meta[match(model_samples, sample_id)]
if (anyNA(model_meta$sample_id) || anyNA(model_meta$inferred_sex) ||
    any(!nzchar(trimws(as.character(model_meta$inferred_sex))))) {
  stop("Model metadata lacks exact inferred-sex coverage", call. = FALSE)
}
model_info <- data.frame(
  group_binary = factor(dge$samples$group_binary[dge$samples$dataset %in% mega], levels = c("Control", "Disease")),
  dataset = factor(dge$samples$dataset[dge$samples$dataset %in% mega]),
  inferred_sex = factor(model_meta$inferred_sex)
)
design_expected <- model.matrix(~ dataset + inferred_sex + group_binary, data = model_info)
if (qr(design_expected)$rank != ncol(design_expected) ||
    !"group_binaryDisease" %in% colnames(design_expected) ||
    !identical(weights$sample_id, model_samples) ||
    any(!is.finite(weights$sample_weight)) || any(weights$sample_weight <= 0) ||
    !identical(design_observed$sample_id, model_samples) ||
    !identical(names(design_observed)[-1L], colnames(design_expected))) {
  stop("Model sample/weight/design identity or estimability failed", call. = FALSE)
}
design_matrix_observed <- as.matrix(
  design_observed[, setdiff(names(design_observed), "sample_id"), with = FALSE]
)
storage.mode(design_matrix_observed) <- "double"
storage.mode(design_expected) <- "double"
dimnames(design_matrix_observed) <- NULL
dimnames(design_expected) <- NULL
# model.matrix() attaches "assign" and "contrasts" bookkeeping that the emitted
# TSV cannot carry. Without stripping them identical() fails on attributes even
# when every cell matches exactly (observed 2026-08-09: 0 of 5,922 cells
# differed, max abs diff 0). This check is about VALUES, so compare values.
attr(design_expected, "assign") <- NULL
attr(design_expected, "contrasts") <- NULL
if (!identical(design_matrix_observed, design_expected)) {
  stop("Emitted model design values differ from the exact reconstructed design", call. = FALSE)
}
expected_model_paths <- list(
  dge_path = normalizePath(file.path(arm_root, "results/integration/merged_dge.rds"), mustWork = TRUE),
  meta_path = normalizePath(file.path(arm_root, "results/integration/meta_matched.rds"), mustWork = TRUE),
  config_path = normalizePath(file.path(project_root, "config/human_datasets.yaml"), mustWork = TRUE),
  gene_metadata_path = normalizePath(file.path(project_root, "data/gencode_v49_gene_metadata.tsv.gz"), mustWork = TRUE)
)
if (nrow(model_manifest) != 1L || model_manifest$n_genes != nrow(dge) ||
    model_manifest$n_samples != length(model_samples) ||
    model_manifest$n_cohorts != length(unique(model_info$dataset)) ||
    model_manifest$n_control != sum(model_info$group_binary == "Control") ||
    model_manifest$n_disease != sum(model_info$group_binary == "Disease") ||
    model_manifest$treat_lfc != 0.25 ||
    any(vapply(names(expected_model_paths), function(field) {
      !identical(model_manifest[[field]][[1L]], expected_model_paths[[field]])
    }, logical(1)))) {
  stop("Model input manifest differs from fitted arm dimensions/TREAT floor", call. = FALSE)
}

# Counts mode must be established by the exact dataset->path->hash mapping, not
# merely by the self-reported preprocessing mode.
override_path <- file.path(
  run_root, "manifests",
  if (identical(mode$counts, "read")) "read_count_overrides.tsv" else "count_overrides.tsv"
)
expected_sources <- fread(override_path)
required_source_fields <- c("dataset", "count_path", "size_bytes", "sha256")
if (!all(required_source_fields %in% names(effective_sources)) ||
    !identical(names(expected_sources), c("dataset", "count_path")) ||
    anyDuplicated(effective_sources$dataset) || anyDuplicated(expected_sources$dataset) ||
    !setequal(effective_sources$dataset, expected_sources$dataset)) {
  stop("Effective count-source manifest has an invalid dataset/schema contract", call. = FALSE)
}
expected_sources <- expected_sources[match(effective_sources$dataset, dataset)]
for (index in seq_len(nrow(effective_sources))) {
  observed_path <- normalizePath(effective_sources$count_path[[index]], mustWork = TRUE)
  expected_path <- normalizePath(expected_sources$count_path[[index]], mustWork = TRUE)
  if (!identical(observed_path, expected_path) ||
      file.info(observed_path)$isdir || nzchar(Sys.readlink(observed_path)) ||
      effective_sources$size_bytes[[index]] != file.info(observed_path)$size ||
      effective_sources$sha256[[index]] != sha256(observed_path)) {
    stop("Effective count source differs from override/path/hash contract: ",
         effective_sources$dataset[[index]], call. = FALSE)
  }
}

curve <- readRDS(file.path(arm_root, "results/integration/voom_curve.rds"))
for (name in c("voom_xy", "voom_line")) {
  value <- curve[[name]]
  if (is.null(value) || is.null(value$x) || is.null(value$y) ||
      length(value$x) == 0L || length(value$x) != length(value$y) ||
      any(!is.finite(value$x)) || any(!is.finite(value$y))) {
    stop("Missing/invalid emitted voom curve: ", name, call. = FALSE)
  }
}

if (arm %in% c("R0", "F_locked")) {
  locked <- readRDS(file.path(run_root, "frozen_sets/locked_merged_dge.rds"))
  if (!identical(colnames(dge), colnames(locked)) ||
      (arm == "F_locked" &&
       (nrow(dge) != 27638L || !identical(rownames(dge), rownames(locked))))) {
    stop(arm, " does not preserve the exact required frozen DGE sample/gene order", call. = FALSE)
  }
}

if (arm == "F_five") {
  legacy_dge <- readRDS(file.path(run_root, "arms/F_legacy/results/integration/merged_dge.rds"))
  legacy_model_samples <- colnames(legacy_dge)[legacy_dge$samples$dataset %in% mega]
  if (!identical(model_samples, legacy_model_samples)) {
    stop("F_five fitted model sample order differs from the F_legacy canonical-five subset", call. = FALSE)
  }
  reused <- c(
    "qc/sample_qc_report.csv",
    "results/integration/meta_matched.rds",
    "results/integration/merged_counts_raw.rds",
    "provenance/effective_count_sources.tsv"
  )
  for (relative in reused) {
    if (sha256(file.path(arm_root, relative)) !=
        sha256(file.path(run_root, "arms/F_legacy", relative))) {
      stop("F_five reuse artifact differs byte-for-byte from F_legacy: ", relative, call. = FALSE)
    }
  }
}

output <- file.path(arm_root, "provenance/arm_validation.tsv")
# NA-safe symlink test: Sys.readlink() yields NA for an absent path and
# nzchar(NA) is TRUE by default, which made this refuse every fresh arm.
if (file.exists(output) ||
    isTRUE(nzchar(Sys.readlink(output), keepNA = TRUE))) {
  stop("Refusing existing arm-validation output", call. = FALSE)
}
fwrite(data.table(
  arm = arm,
  status = "PASS",
  raw_genes = nrow(raw),
  dge_genes = nrow(dge),
  dge_samples = ncol(dge),
  model_samples = length(model_samples),
  treat_lfc = 0.25
), output, sep = "\t")
cat("PASS exact per-arm preprocessing/model/numerical contract: ", arm, "\n", sep = "")
