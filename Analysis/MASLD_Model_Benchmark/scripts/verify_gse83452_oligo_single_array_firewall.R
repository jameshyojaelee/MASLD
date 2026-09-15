#!/usr/bin/env Rscript

# Prove on real GSE83452 CEL files that no held array's summary can see another
# held array, on the GPL16686 code path.
#
# This is deliberately a separate script from the GPL570 fixture. The GPL570
# pass is evidence about affy/frma/rwaFit2 and says nothing about oligo or
# SCAN.UPC, which are different packages with different summarization code. The
# hazard is assumed present here until measured absent.
#
# The prior is that it IS present and worse: oligo::rma summarizes by median
# polish over the probe-by-array matrix, which borrows across arrays by
# construction, where fRMA at least anchors on frozen probe effects.
#
# Two families of context, exactly as on GPL570:
#
#   pipeline contexts  One CEL per summarization call, differing only in what
#                      else the run processed and in what order. These MUST be
#                      bit-identical. This is the separation.
#
#   estimator contexts The same target inside a multi-array call. Expected to
#                      diverge; recorded as the standing hazard measurement.
#
# The 1e-12 perturbation control and the distinct-target-digest check are
# carried over unchanged, so a green result cannot be a broken comparator or a
# collapsed summarizer.
#
# SCAN is probed separately. The runtime capability probe recorded it as failing
# with 'could not find function "%dopar%"'. foreach is installed in the frozen
# runtime, so that is an attach problem in SCAN.UPC rather than a missing
# dependency; a sequential backend is registered and the result recorded either
# way. Nothing is installed.
#
# Targets and companions are chosen by position in the sorted GSM accession
# list. No label, phenotype, or intensity statistic takes part. This script
# reads no label, reads no GEO series matrix, runs no all-sample RMA, runs no
# across-array normalization, exports no expression value, and fits no model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: verify_gse83452_oligo_single_array_firewall.R CEL_DIR MANIFEST_TSV OUTPUT_DIR")
}
cel_dir <- args[[1L]]
manifest_path <- args[[2L]]
output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages({
  library(oligo)
  library(pd.hugene.2.0.st)
  library(foreach)
  library(jsonlite)
})
# SCAN.UPC uses %dopar% without importing it. Register a sequential backend so
# the probe measures SCAN itself rather than a packaging defect.
foreach::registerDoSEQ()

EXPECTED <- c(
  "oligo" = "1.74.0", "affyio" = "1.80.0", "affxparser" = "1.82.0",
  "pd.hugene.2.0.st" = "3.14.1", "SCAN.UPC" = "2.52.0"
)
observed <- vapply(
  names(EXPECTED), function(package) as.character(packageVersion(package)), character(1)
)
if (!identical(unname(observed), unname(EXPECTED))) {
  stop(sprintf(
    "GPL16686 summarization package versions differ: %s",
    paste(sprintf("%s=%s", names(observed), observed), collapse = ",")
  ))
}

## -------------------------------------------------------- frozen CEL manifest
manifest <- utils::read.delim(
  manifest_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
manifest <- manifest[manifest$series == "GSE83452", , drop = FALSE]
manifest <- manifest[order(manifest$sample_accession), , drop = FALSE]
if (nrow(manifest) != 231L) stop("GSE83452 CEL manifest is not 231 records")
accessions <- manifest$sample_accession
cel_path <- function(accession) file.path(cel_dir, paste0(accession, ".CEL"))
if (!all(file.exists(cel_path(accessions)))) stop("a staged GSE83452 CEL is absent")

# Label-blind, position-based selection. Companion sets are smaller than on
# GPL570 because HuGene 2.0 ST arrays carry 2.6M probes each.
target_positions <- c(1L, 116L, 231L)
companion_p_positions <- 2L:4L
companion_q_positions <- 228L:230L
targets <- accessions[target_positions]
companion_p <- accessions[companion_p_positions]
companion_q <- accessions[companion_q_positions]
if (length(intersect(targets, c(companion_p, companion_q))) != 0L) {
  stop("target and companion arrays overlap")
}

digest_doubles <- function(values) {
  path <- tempfile(fileext = ".f8")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, "wb")
  writeBin(as.double(values), connection, size = 8L, endian = "little")
  close(connection)
  substr(system2("sha256sum", shQuote(path), stdout = TRUE), 1L, 64L)
}

# The feature axis is whatever the core summarizer returns; it is captured on
# the first call and every later call must match it exactly.
feature_axis <- NULL

# THE production call: exactly one CEL, normalization off, core target.
summarize_one <- function(accession) {
  feature_set <- oligo::read.celfiles(cel_path(accession), verbose = FALSE)
  summary <- oligo::rma(feature_set, background = TRUE, normalize = FALSE, target = "core")
  values <- Biobase::exprs(summary)
  if (ncol(values) != 1L) stop(sprintf("%s was not summarized as a single array", accession))
  if (is.null(feature_axis)) {
    feature_axis <<- rownames(values)
  } else if (!identical(rownames(values), feature_axis)) {
    stop(sprintf("%s core feature axis differs from the first array", accession))
  }
  result <- as.numeric(values[feature_axis, 1L])
  rm(feature_set, summary, values)
  invisible(gc(verbose = FALSE))
  result
}

# A multi-array call. Retained only to keep measuring the hazard; never a path.
summarize_in_batch <- function(target, accession_set) {
  feature_set <- oligo::read.celfiles(cel_path(accession_set), verbose = FALSE)
  summary <- oligo::rma(feature_set, background = TRUE, normalize = FALSE, target = "core")
  values <- Biobase::exprs(summary)
  if (ncol(values) != length(accession_set)) {
    stop("oligo column axis differs from the requested array set")
  }
  if (!identical(rownames(values), feature_axis)) {
    stop("multi-array core feature axis differs from the single-array axis")
  }
  result <- as.numeric(values[feature_axis, match(target, accession_set)])
  rm(feature_set, summary, values)
  invisible(gc(verbose = FALSE))
  result
}

pipeline_contexts <- list(
  solo_only = function(target) summarize_one(target),
  solo_after_P = function(target) {
    for (companion in companion_p) invisible(summarize_one(companion))
    summarize_one(target)
  },
  solo_before_P = function(target) {
    value <- summarize_one(target)
    for (companion in companion_p) invisible(summarize_one(companion))
    value
  },
  solo_interleaved_Q = function(target) {
    invisible(summarize_one(companion_q[[1L]]))
    value <- summarize_one(target)
    for (companion in companion_q[2L:length(companion_q)]) invisible(summarize_one(companion))
    value
  }
)
estimator_contexts <- list(
  batch_P_first = function(target) summarize_in_batch(target, c(target, companion_p)),
  batch_P_last = function(target) summarize_in_batch(target, c(companion_p, target)),
  batch_Q_first = function(target) summarize_in_batch(target, c(target, companion_q))
)

compare_to <- function(reference, values) {
  list(
    bitwise_identical = identical(values, reference),
    differing_elements = sum(values != reference),
    maximum_absolute_difference = max(abs(values - reference))
  )
}

target_records <- list()
alone_digests <- character(0)
for (target in targets) {
  pipeline <- lapply(pipeline_contexts, function(context) context(target))
  reference <- pipeline[["solo_only"]]
  if (!all(is.finite(reference))) stop("single-array summary carries a nonfinite value")

  pipeline_digests <- vapply(pipeline, digest_doubles, character(1))
  pipeline_comparison <- lapply(pipeline, function(values) compare_to(reference, values))

  estimator <- lapply(estimator_contexts, function(context) context(target))
  estimator_digests <- vapply(estimator, digest_doubles, character(1))
  estimator_comparison <- lapply(estimator, function(values) compare_to(reference, values))

  perturbed <- reference
  perturbed[[1L]] <- perturbed[[1L]] + 1e-12
  perturbed_digest <- digest_doubles(perturbed)

  quantiles <- unname(stats::quantile(reference, probs = c(0, 0.25, 0.5, 0.75, 1)))
  target_records[[target]] <- list(
    sample_accession = target,
    pipeline_contexts_compared = length(pipeline),
    pipeline_digests = as.list(pipeline_digests),
    pipeline_all_digests_identical = length(unique(unname(pipeline_digests))) == 1L,
    pipeline_all_bitwise_identical =
      all(vapply(pipeline_comparison, function(x) isTRUE(x$bitwise_identical), logical(1))),
    pipeline_differing_elements = lapply(pipeline_comparison, function(x) x$differing_elements),
    pipeline_maximum_absolute_difference =
      lapply(pipeline_comparison, function(x) x$maximum_absolute_difference),
    pipeline_arrays_per_summarization_call = 1L,
    estimator_contexts_compared = length(estimator),
    estimator_digests = as.list(estimator_digests),
    estimator_all_bitwise_identical =
      all(vapply(estimator_comparison, function(x) isTRUE(x$bitwise_identical), logical(1))),
    estimator_differing_elements = lapply(estimator_comparison, function(x) x$differing_elements),
    estimator_maximum_absolute_difference =
      lapply(estimator_comparison, function(x) x$maximum_absolute_difference),
    perturbation_control_digest_changed =
      !identical(perturbed_digest, pipeline_digests[["solo_only"]]),
    perturbation_control_magnitude = 1e-12,
    features = length(reference),
    finite_values = sum(is.finite(reference)),
    min = quantiles[[1L]], q25 = quantiles[[2L]], median = quantiles[[3L]],
    q75 = quantiles[[4L]], max = quantiles[[5L]]
  )
  alone_digests <- c(alone_digests, unname(pipeline_digests[["solo_only"]]))
  rm(pipeline, estimator)
  invisible(gc(verbose = FALSE))
}

## ------------------------------------------------- SCAN capability and hazard
scan_probe <- local({
  attempt <- function(expression) {
    tryCatch(
      list(status = "pass", error = "not_applicable", value = expression),
      error = function(condition) list(status = "failed", error = conditionMessage(condition), value = NULL)
    )
  }
  target <- targets[[1L]]
  single <- attempt({
    values <- Biobase::exprs(SCAN.UPC::SCAN(cel_path(target), outFilePath = NA, verbose = FALSE))
    as.numeric(values[, 1L])
  })
  record <- list(
    method = "SCAN_with_pd.hugene.2.0.st",
    foreach_sequential_backend_registered = TRUE,
    single_array_status = single$status,
    single_array_error = single$error
  )
  if (identical(single$status, "pass")) {
    record$single_array_features <- length(single$value)
    record$single_array_finite <- sum(is.finite(single$value))
    repeated <- attempt({
      values <- Biobase::exprs(SCAN.UPC::SCAN(cel_path(target), outFilePath = NA, verbose = FALSE))
      as.numeric(values[, 1L])
    })
    if (identical(repeated$status, "pass")) {
      record$repeat_call_bitwise_identical <- identical(repeated$value, single$value)
    }
    batched <- attempt({
      values <- Biobase::exprs(SCAN.UPC::SCAN(
        cel_path(c(target, companion_p)), outFilePath = NA, verbose = FALSE
      ))
      as.numeric(values[, 1L])
    })
    record$batched_status <- batched$status
    record$batched_error <- batched$error
    if (identical(batched$status, "pass")) {
      record$batched_vs_single_bitwise_identical <- identical(batched$value, single$value)
      record$batched_vs_single_differing_elements <- sum(batched$value != single$value)
      record$batched_vs_single_maximum_absolute_difference <-
        max(abs(batched$value - single$value))
    }
  }
  record
})

pipeline_pass <- all(vapply(
  target_records, function(record) isTRUE(record$pipeline_all_bitwise_identical), logical(1)
))
receipt <- list(
  schema_version = "masld-bench-gse83452-oligo-single-array-firewall-v1",
  status = if (pipeline_pass) "pass_single_array_isolation_on_real_CEL" else "fail_single_array_isolation_on_real_CEL",
  series = "GSE83452",
  platform_id = "GPL16686",
  cohort_family_id = "antwerp_inserm_shared",
  method = "single_array_core_summary_normalize_false",
  summarize_rule = "oligo_rma_median_polish_core_target",
  normalize_rule = "none",
  background_rule = "rma",
  held_array_application = "one_array_at_a_time_with_frozen_training_object",
  package_versions = as.list(observed),
  core_features = if (is.null(feature_axis)) NA_integer_ else length(feature_axis),
  cel_records_available = nrow(manifest),
  target_selection = "sorted_GSM_accession_positions_1_116_231",
  companion_selection = "sorted_GSM_accession_positions_2_to_4_and_228_to_230",
  selection_used_labels = FALSE,
  selection_used_intensity = FALSE,
  targets = as.list(targets),
  companion_set_P = as.list(companion_p),
  companion_set_Q = as.list(companion_q),
  pipeline_context_names = as.list(names(pipeline_contexts)),
  estimator_context_names = as.list(names(estimator_contexts)),
  per_target = unname(target_records),
  every_target_pipeline_bitwise_identical = pipeline_pass,
  every_target_pipeline_digest_identical = all(vapply(
    target_records, function(record) isTRUE(record$pipeline_all_digests_identical), logical(1)
  )),
  every_perturbation_control_detected = all(vapply(
    target_records, function(record) isTRUE(record$perturbation_control_digest_changed), logical(1)
  )),
  distinct_targets_have_distinct_digests = length(unique(alone_digests)) == length(alone_digests),
  multi_array_call_is_batch_dependent = !all(vapply(
    target_records, function(record) isTRUE(record$estimator_all_bitwise_identical), logical(1)
  )),
  multi_array_call_used_as_preprocessing_path = FALSE,
  gpl570_evidence_carried_over = FALSE,
  scan_probe = scan_probe,
  all_sample_RMA_run = FALSE,
  across_array_quantile_normalization_run = FALSE,
  GEO_series_matrix_read = FALSE,
  expression_values_exported = FALSE,
  labels_read = FALSE,
  model_training_activated = FALSE,
  sealed_outcomes_read = FALSE
)
writeLines(
  jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE, digits = NA),
  file.path(output_dir, "oligo_single_array_firewall.json")
)
writeLines(
  capture.output(print(utils::sessionInfo())),
  file.path(output_dir, "R_sessionInfo.txt")
)
cat(sprintf(
  "targets=%d core_features=%s pipeline_bitwise_identical=%s perturbation_detected=%s multi_array_batch_dependent=%s scan_single=%s\n",
  length(targets),
  if (is.null(feature_axis)) "NA" else length(feature_axis),
  receipt$every_target_pipeline_bitwise_identical,
  receipt$every_perturbation_control_detected,
  receipt$multi_array_call_is_batch_dependent,
  scan_probe$single_array_status
))
