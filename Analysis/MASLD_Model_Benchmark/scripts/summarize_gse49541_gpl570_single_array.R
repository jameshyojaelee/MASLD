#!/usr/bin/env Rscript

# Label-blind single-array fRMA summarization of all 72 GSE49541 GPL570 arrays.
#
# Every array is read and summarized on its own. The AffyBatch never holds more
# than one array, so the frozen hgu133plus2frmavecs normalization target and
# probe effects are the only cross-array information that ever touches a value,
# and those come from the vector package, not from this cohort. That property is
# proven separately, on these same CEL files, by
# scripts/verify_gse49541_frma_single_array_firewall.R.
#
# Two matrices leave this script. The platform-feature matrix is the frozen
# 54675-row GPL570 axis. The gene matrix keeps only features with one unambiguous
# Entrez ID and one unambiguous GENCODE v49 Ensembl gene, collapsed by the
# training-frozen median of the log-scale features. The collapse is a within-array
# median over a frozen feature set, so it is per-array independent as well.
#
# QC is descriptive. Per-array metrics are computed from that array alone.
# Cross-array metrics (RLE, PCA distance) are computed only after both matrices
# are final, are recorded as descriptive flags, never modify a value, and never
# remove a sample: the requirements set automatic_exclusion_from_descriptive_flag
# to false and label_blind to true.
#
# This script reads no label, reads no GEO series matrix, runs no all-sample RMA,
# and runs no across-array normalization.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: summarize_gse49541_gpl570_single_array.R CEL_DIR MANIFEST_TSV AXIS_TSV CROSSWALK_TSV OUTPUT_DIR")
}
cel_dir <- args[[1L]]
manifest_path <- args[[2L]]
axis_path <- args[[3L]]
crosswalk_path <- args[[4L]]
output_dir <- args[[5L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages({
  library(affy)
  library(affyio)
  library(frma)
  library(jsonlite)
})

EXPECTED <- c(
  "affy" = "1.88.0", "affyio" = "1.80.0", "frma" = "1.62.0",
  "hgu133plus2cdf" = "2.18.0", "hgu133plus2frmavecs" = "1.5.0"
)
observed <- vapply(
  names(EXPECTED), function(package) as.character(packageVersion(package)), character(1)
)
if (!identical(unname(observed), unname(EXPECTED))) {
  stop(sprintf(
    "GPL570 summarization package versions differ: %s",
    paste(sprintf("%s=%s", names(observed), observed), collapse = ",")
  ))
}

sha256_of <- function(path) {
  substr(system2("sha256sum", shQuote(path), stdout = TRUE), 1L, 64L)
}

## ------------------------------------------------------- frozen platform axis
axis <- utils::read.delim(
  axis_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
platform_features <- axis[["ID"]]
if (length(platform_features) != 54675L || anyDuplicated(platform_features) > 0L) {
  stop("frozen GPL570 platform axis is not 54675 unique features")
}

## -------------------------------------------------------- frozen CEL manifest
manifest <- utils::read.delim(
  manifest_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
manifest <- manifest[order(manifest$sample_accession), , drop = FALSE]
if (nrow(manifest) != 72L) stop("GSE49541 CEL manifest is not 72 records")
accessions <- manifest$sample_accession
cel_path <- function(accession) file.path(cel_dir, paste0(accession, ".CEL"))
if (!all(file.exists(cel_path(accessions)))) stop("a staged GSE49541 CEL is absent")

## ------------------------------------------------- frozen GENCODE v49 crosswalk
crosswalk <- utils::read.delim(
  crosswalk_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
if (nrow(crosswalk) != 54675L || !identical(crosswalk$platform_feature_id, platform_features)) {
  stop("frozen crosswalk axis differs from the frozen GPL570 platform axis")
}
eligible <- crosswalk[crosswalk$eligible_for_gene_matrix == "true", , drop = FALSE]
if (nrow(eligible) < 1L) stop("frozen crosswalk retains no one-to-one feature")
gene_ids <- sort(unique(eligible$ensembl_gene_id))
feature_row <- match(eligible$platform_feature_id, platform_features)
gene_index <- match(eligible$ensembl_gene_id, gene_ids)
gene_groups <- split(feature_row, gene_index)
names(gene_groups) <- gene_ids[as.integer(names(gene_groups))]
gene_groups <- gene_groups[gene_ids]

## ------------------------------------------------- frozen fRMA vector package
loader <- new.env(parent = emptyenv())
utils::data("hgu133plus2frmavecs", package = "hgu133plus2frmavecs", envir = loader)
frma_vectors <- loader[["hgu133plus2frmavecs"]]
vector_digest <- local({
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(frma_vectors, path, compress = FALSE)
  sha256_of(path)
})

## ------------------------------------------ one array at a time, never a batch
platform_matrix <- matrix(
  NA_real_, nrow = length(platform_features), ncol = length(accessions),
  dimnames = list(platform_features, accessions)
)
gene_matrix <- matrix(
  NA_real_, nrow = length(gene_ids), ncol = length(accessions),
  dimnames = list(gene_ids, accessions)
)
per_array <- vector("list", length(accessions))

for (index in seq_along(accessions)) {
  accession <- accessions[[index]]
  path <- cel_path(accession)

  header <- affyio::read.celfile.header(path, info = "full")
  dimensions <- header[["CEL dimensions"]]
  scan_date <- if ("ScanDate" %in% names(header)) as.character(header[["ScanDate"]]) else "not_available"

  batch <- affy::ReadAffy(filenames = path)
  if (length(Biobase::sampleNames(batch)) != 1L) {
    stop(sprintf("%s was not read as a single array", accession))
  }
  raw <- as.numeric(affy::intensity(batch))
  raw_quantiles <- unname(stats::quantile(raw, probs = c(0, 0.25, 0.5, 0.75, 1)))

  degradation <- tryCatch(
    {
      fit <- affy::AffyRNAdeg(batch)
      list(slope = as.numeric(fit$slope)[[1L]], pvalue = as.numeric(fit$pvalue)[[1L]])
    },
    error = function(condition) list(slope = NA_real_, pvalue = NA_real_)
  )

  summary <- frma::frma(batch, input.vecs = frma_vectors, verbose = FALSE)
  values <- Biobase::exprs(summary)
  if (!setequal(rownames(values), platform_features) || ncol(values) != 1L) {
    stop(sprintf("%s fRMA summary axis differs from the frozen GPL570 axis", accession))
  }
  column <- as.numeric(values[platform_features, 1L])
  if (!all(is.finite(column))) stop(sprintf("%s summary carries a nonfinite value", accession))
  platform_matrix[, index] <- column

  # Within-array median collapse over the frozen one-to-one feature set.
  gene_matrix[, index] <- vapply(
    gene_groups, function(rows) stats::median(column[rows]), numeric(1)
  )

  summary_quantiles <- unname(stats::quantile(column, probs = c(0, 0.25, 0.5, 0.75, 1)))
  per_array[[index]] <- list(
    sample_accession = accession,
    arrays_read_in_this_call = 1L,
    cdf_name = as.character(header[["cdfName"]]),
    grid_columns = as.integer(dimensions[[1L]]),
    grid_rows = as.integer(dimensions[[2L]]),
    scan_date = scan_date,
    cel_parse_ok = TRUE,
    array_type_matches_GPL570 = identical(as.character(header[["cdfName"]]), "HG-U133_Plus_2"),
    raw_probes = length(raw),
    raw_nonpositive = sum(raw <= 0),
    raw_min = raw_quantiles[[1L]], raw_q25 = raw_quantiles[[2L]],
    raw_median = raw_quantiles[[3L]], raw_q75 = raw_quantiles[[4L]],
    raw_max = raw_quantiles[[5L]],
    rna_degradation_slope = degradation$slope,
    rna_degradation_pvalue = degradation$pvalue,
    summary_features = length(column),
    summary_finite = sum(is.finite(column)),
    summary_min = summary_quantiles[[1L]], summary_q25 = summary_quantiles[[2L]],
    summary_median = summary_quantiles[[3L]], summary_q75 = summary_quantiles[[4L]],
    summary_max = summary_quantiles[[5L]]
  )
  rm(batch, summary, values, raw)
  invisible(gc(verbose = FALSE))
}

if (any(!is.finite(platform_matrix))) stop("platform feature matrix carries a nonfinite value")
if (any(!is.finite(gene_matrix))) stop("gene matrix carries a nonfinite value")

## ------------------------- descriptive cross-array QC, computed after the fact
# These never modify a value and never remove a sample. The requirements set
# automatic_exclusion_from_descriptive_flag = false, and
# held_cohort_intensity_distribution_visible_to_training = false, so these
# numbers are for human audit only and may not enter a preprocessing object.
feature_median <- apply(platform_matrix, 1L, stats::median)
rle <- platform_matrix - feature_median
rle_median <- apply(rle, 2L, stats::median)
rle_iqr <- apply(rle, 2L, stats::IQR)

centered <- t(gene_matrix - rowMeans(gene_matrix))
decomposition <- stats::prcomp(centered, center = FALSE, scale. = FALSE)
components <- decomposition$x[, seq_len(min(5L, ncol(decomposition$x))), drop = FALSE]
component_centers <- apply(components, 2L, stats::median)
component_scales <- apply(components, 2L, function(values) {
  scale <- stats::mad(values)
  if (!is.finite(scale) || scale <= 0) 1 else scale
})
robust_distance <- sqrt(rowSums(
  (sweep(sweep(components, 2L, component_centers, "-"), 2L, component_scales, "/"))^2
))
variance_explained <- (decomposition$sdev^2) / sum(decomposition$sdev^2)

for (index in seq_along(accessions)) {
  per_array[[index]]$rle_median <- unname(rle_median[[index]])
  per_array[[index]]$rle_iqr <- unname(rle_iqr[[index]])
  per_array[[index]]$pca_robust_distance <- unname(robust_distance[[index]])
  per_array[[index]]$descriptive_flag_only <- TRUE
  per_array[[index]]$excluded <- FALSE
}

## ------------------------------------------------------------------- outputs
# The whole lane now rests on bit-exactness, so the written matrix is read back
# and required to be bitwise identical to what is in memory. That turns the
# usual assumption about text round-tripping of IEEE-754 doubles into a
# measurement, and fails the job rather than the analysis if it does not hold.
write_matrix <- function(values, id_column, path) {
  # write.table converts doubles at 15 significant digits, which is short of the
  # 17 an IEEE-754 double needs to survive a text round trip; on a log2
  # expression scale roughly 93 percent of values lose their last bits. The
  # values are therefore formatted at 17 significant digits explicitly. %g still
  # drops trailing zeros, so exactly representable values stay short.
  text <- matrix(
    sprintf("%.17g", values), nrow = nrow(values), ncol = ncol(values)
  )
  frame <- data.frame(
    id = rownames(values), text, check.names = FALSE, stringsAsFactors = FALSE
  )
  names(frame) <- c(id_column, colnames(values))
  utils::write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  restored <- utils::read.delim(
    path, sep = "\t", quote = "", check.names = FALSE, row.names = 1L
  )
  restored <- as.matrix(restored)
  if (!identical(dim(restored), dim(values)) ||
      !identical(rownames(restored), rownames(values)) ||
      !identical(colnames(restored), colnames(values)) ||
      !identical(as.numeric(restored), as.numeric(values))) {
    stop(sprintf("written matrix does not round-trip bitwise: %s", path))
  }
  invisible(TRUE)
}
platform_path <- file.path(output_dir, "gse49541_gpl570_platform_feature_matrix.tsv")
gene_path <- file.path(output_dir, "gse49541_gpl570_gene_matrix.tsv")
write_matrix(platform_matrix, "platform_feature_id", platform_path)
write_matrix(gene_matrix, "ensembl_gene_id", gene_path)
writeLines(accessions, file.path(output_dir, "summarized_columns.txt"))

utils::write.table(
  do.call(rbind, lapply(per_array, function(record) as.data.frame(record, stringsAsFactors = FALSE))),
  file.path(output_dir, "gse49541_qc_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

# The preprocessing object is exactly what a held array is transformed by: the
# frozen fRMA vectors, the frozen crosswalk, and the frozen collapse rule.
preprocessing_object <- list(
  method = "frma_with_exact_hgu133plus2frmavecs",
  summarize_rule = "robust_weighted_average",
  normalize_rule = "quantile_to_frozen_hgu133plus2frmavecs_target",
  background_rule = "rma",
  target = "probeset",
  frma_vector_object_sha256 = vector_digest,
  crosswalk_sha256 = sha256_of(crosswalk_path),
  platform_axis_sha256 = sha256_of(axis_path),
  gene_collapse_rule = "training_frozen_median_of_log_scale_features_after_one_to_one_mapping",
  package_versions = as.list(observed)
)
preprocessing_object_path <- file.path(output_dir, "preprocessing_object.json")
writeLines(
  jsonlite::toJSON(preprocessing_object, auto_unbox = TRUE, pretty = TRUE, digits = NA),
  preprocessing_object_path
)

receipt <- list(
  schema_version = "masld-bench-gse49541-gpl570-single-array-summary-v1",
  status = "pass_label_blind_single_array_summarization",
  series = "GSE49541",
  platform_id = "GPL570",
  cohort_family_id = "gse31803_gse49541_fibrosis_array",
  method = "frma_with_exact_hgu133plus2frmavecs",
  arrays_read_per_frma_call = 1L,
  arrays_summarized = length(accessions),
  platform_features = nrow(platform_matrix),
  gene_features = nrow(gene_matrix),
  eligible_platform_features = nrow(eligible),
  package_versions = as.list(observed),
  preprocessing_object = preprocessing_object,
  preprocessing_object_sha256 = sha256_of(preprocessing_object_path),
  matrices_round_trip_bitwise = TRUE,
  platform_feature_matrix_sha256 = sha256_of(platform_path),
  gene_matrix_sha256 = sha256_of(gene_path),
  hard_failures = list(
    CEL_parse_failure = 0L,
    array_type_mismatch = sum(!vapply(per_array, function(r) isTRUE(r$array_type_matches_GPL570), logical(1))),
    sample_axis_mismatch = 0L,
    nonfinite_summary = 0L,
    duplicate_feature_or_sample_ID = as.integer(
      anyDuplicated(rownames(platform_matrix)) > 0L || anyDuplicated(colnames(platform_matrix)) > 0L
    ),
    no_positive_intensity = sum(vapply(per_array, function(r) r$raw_nonpositive == r$raw_probes, logical(1)))
  ),
  pca_variance_explained = as.list(unname(variance_explained[seq_len(min(5L, length(variance_explained)))])),
  descriptive_flags_computed_after_matrices_were_final = TRUE,
  automatic_exclusion_from_descriptive_flag = FALSE,
  cross_array_QC_visible_to_training = FALSE,
  label_blind = TRUE,
  samples_excluded = 0L,
  per_array = per_array,
  all_sample_RMA_run = FALSE,
  across_array_quantile_normalization_run = FALSE,
  GEO_series_matrix_read = FALSE,
  labels_read = FALSE,
  model_training_activated = FALSE,
  sealed_outcomes_read = FALSE
)
writeLines(
  jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE, digits = NA),
  file.path(output_dir, "summarization_receipt.json")
)
writeLines(
  capture.output(print(utils::sessionInfo())),
  file.path(output_dir, "R_sessionInfo.txt")
)
utils::write.table(
  data.frame(package = names(observed), version = unname(observed), stringsAsFactors = FALSE),
  file.path(output_dir, "package_versions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
cat(sprintf(
  "arrays=%d platform_features=%d genes=%d\n",
  length(accessions), nrow(platform_matrix), nrow(gene_matrix)
))
