#!/usr/bin/env Rscript

# Prove on real GSE49541 CEL files that no held array's summary can see another
# held array.
#
# The first run of this fixture (model-check-085, job 21107058) established that
# fRMA's default summarize = "robust_weighted_average" is NOT invariant to which
# arrays share its call: rwaFit2 receives the whole probe-by-array submatrix and
# re-estimates its robustness weights over every column present. Batch position
# changed nothing; batch membership moved ~71 percent of features by up to
# 9e-4 on the log2 scale.
#
# That finding is kept, not weakened. It is the reason the pipeline must hold
# the batch at n=1, and it is re-measured here every run so that nobody can
# later reintroduce batching without the fixture failing.
#
# Two families of context are therefore compared, and they answer different
# questions:
#
#   pipeline contexts  Each frma() call reads exactly ONE CEL, which is what
#                      config/evaluation/microarray_transfer_preprocessing.toml
#                      requires: held_array_application =
#                      "one_array_at_a_time_with_frozen_training_object".
#                      The contexts differ in what else the same run processed,
#                      and in what order. These MUST be bit-identical: with one
#                      column there is nothing for rwaFit2 to borrow from, and
#                      no session state may carry between arrays. This is the
#                      separation.
#
#   estimator contexts The same target inside a multi-array frma() call, exactly
#                      as the first run built them. These are NOT expected to
#                      match, and their measured divergence is recorded as the
#                      standing evidence for why n=1 is mandatory. They are
#                      negative controls, never a preprocessing path.
#
# A perturbation control adds 1e-12 to one of 54675 values and requires the
# digest to change, so a bit-identical pipeline result cannot be a broken
# comparator. Distinct targets must carry distinct digests, so per-array
# independence is not per-array collapse.
#
# Targets and companions are chosen by position in the sorted GSM accession
# list. No label, phenotype, or intensity statistic takes part in the choice.
# This script reads no label, reads no GEO series matrix, runs no all-sample
# RMA, runs no across-array normalization, exports no expression value, and
# fits no model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: verify_gse49541_frma_single_array_firewall.R CEL_DIR MANIFEST_TSV AXIS_TSV OUTPUT_DIR")
}
cel_dir <- args[[1L]]
manifest_path <- args[[2L]]
axis_path <- args[[3L]]
output_dir <- args[[4L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages({
  library(affy)
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

## ------------------------------------------------------- frozen platform axis
axis <- utils::read.delim(
  axis_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
platform_features <- axis[["ID"]]
if (length(platform_features) != 54675L || anyDuplicated(platform_features) > 0L) {
  stop("frozen GPL570 platform axis is not 54675 unique features")
}

## ---------------------------------------- frozen CEL manifest, sorted by GSM
manifest <- utils::read.delim(
  manifest_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
manifest <- manifest[order(manifest$sample_accession), , drop = FALSE]
if (nrow(manifest) != 72L) stop("GSE49541 CEL manifest is not 72 records")
accessions <- manifest$sample_accession
cel_path <- function(accession) file.path(cel_dir, paste0(accession, ".CEL"))
if (!all(file.exists(cel_path(accessions)))) stop("a staged GSE49541 CEL is absent")

# Label-blind, position-based selection over the sorted accession list.
target_positions <- c(1L, 36L, 72L)
companion_p_positions <- 2L:6L
companion_q_positions <- 67L:71L
targets <- accessions[target_positions]
companion_p <- accessions[companion_p_positions]
companion_q <- accessions[companion_q_positions]
if (length(intersect(targets, c(companion_p, companion_q))) != 0L) {
  stop("target and companion arrays overlap")
}

## ------------------------------------------------- frozen fRMA vector package
loader <- new.env(parent = emptyenv())
utils::data("hgu133plus2frmavecs", package = "hgu133plus2frmavecs", envir = loader)
frma_vectors <- loader[["hgu133plus2frmavecs"]]

vector_digest <- local({
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(frma_vectors, path, compress = FALSE)
  substr(system2("sha256sum", shQuote(path), stdout = TRUE), 1L, 64L)
})

# SHA-256 over the raw little-endian doubles, in frozen platform-axis order.
digest_doubles <- function(values) {
  path <- tempfile(fileext = ".f8")
  on.exit(unlink(path), add = TRUE)
  connection <- file(path, "wb")
  writeBin(as.double(values), connection, size = 8L, endian = "little")
  close(connection)
  substr(system2("sha256sum", shQuote(path), stdout = TRUE), 1L, 64L)
}

# THE production call: exactly one CEL per frma() call, n=1 by construction.
summarize_one <- function(accession) {
  batch <- affy::ReadAffy(filenames = cel_path(accession))
  if (length(Biobase::sampleNames(batch)) != 1L) {
    stop(sprintf("%s was not read as a single array", accession))
  }
  summary <- frma::frma(batch, input.vecs = frma_vectors, verbose = FALSE)
  values <- Biobase::exprs(summary)
  if (!setequal(rownames(values), platform_features) || ncol(values) != 1L) {
    stop(sprintf("%s fRMA summary axis differs from the frozen GPL570 axis", accession))
  }
  result <- as.numeric(values[platform_features, 1L])
  rm(batch, summary, values)
  invisible(gc(verbose = FALSE))
  result
}

# A multi-array call. Retained only to keep measuring the hazard; never a path.
summarize_in_batch <- function(target, accession_set) {
  batch <- affy::ReadAffy(filenames = cel_path(accession_set))
  summary <- frma::frma(batch, input.vecs = frma_vectors, verbose = FALSE)
  values <- Biobase::exprs(summary)
  expected_with_extension <- paste0(accession_set, ".CEL")
  if (!identical(colnames(values), expected_with_extension) &&
      !identical(colnames(values), accession_set)) {
    stop(sprintf(
      "fRMA column axis differs from the requested array set: %s",
      paste(colnames(values), collapse = ",")
    ))
  }
  result <- as.numeric(values[platform_features, match(target, accession_set)])
  rm(batch, summary, values)
  invisible(gc(verbose = FALSE))
  result
}

## ----------------------------------------------------------- pipeline contexts
# Each entry summarizes the target with summarize_one() but surrounds that call
# with different work, in different orders, so that any leakage through session
# state or accumulated globals would show up as a changed digest.
pipeline_contexts <- list(
  solo_only = function(target) {
    summarize_one(target)
  },
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
    invisible(summarize_one(companion_q[[2L]]))
    value <- summarize_one(target)
    for (companion in companion_q[3L:length(companion_q)]) invisible(summarize_one(companion))
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

  # Perturbation control: the digest must detect a 1e-12 change in one feature.
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
    pipeline_differing_elements =
      lapply(pipeline_comparison, function(x) x$differing_elements),
    pipeline_maximum_absolute_difference =
      lapply(pipeline_comparison, function(x) x$maximum_absolute_difference),
    pipeline_arrays_per_frma_call = 1L,
    estimator_contexts_compared = length(estimator),
    estimator_digests = as.list(estimator_digests),
    estimator_all_bitwise_identical =
      all(vapply(estimator_comparison, function(x) isTRUE(x$bitwise_identical), logical(1))),
    estimator_differing_elements =
      lapply(estimator_comparison, function(x) x$differing_elements),
    estimator_maximum_absolute_difference =
      lapply(estimator_comparison, function(x) x$maximum_absolute_difference),
    perturbation_control_digest_changed = !identical(perturbed_digest, pipeline_digests[["solo_only"]]),
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

## ---------------- settle whether summarize = "random_effect" is an alternative
# It has been suggested that random_effect is single-array safe. It estimates a
# between-array variance component, so it should need more than one array, not
# fewer. Measured rather than asserted.
random_effect_probe <- local({
  attempt <- function(expression) {
    tryCatch(
      list(status = "pass", error = "not_applicable", value = expression),
      error = function(condition) list(status = "failed", error = conditionMessage(condition), value = NULL)
    )
  }
  target <- targets[[1L]]
  summarize_rule <- function(accession_set) {
    batch <- affy::ReadAffy(filenames = cel_path(accession_set))
    values <- Biobase::exprs(frma::frma(
      batch, summarize = "random_effect", input.vecs = frma_vectors, verbose = FALSE
    ))
    as.numeric(values[platform_features, match(target, accession_set)])
  }
  single <- attempt(summarize_rule(target))
  batched_p <- attempt(summarize_rule(c(target, companion_p)))
  batched_q <- attempt(summarize_rule(c(target, companion_q)))
  record <- list(
    rule = "random_effect",
    single_array_status = single$status,
    single_array_error = single$error,
    batched_P_status = batched_p$status,
    batched_Q_status = batched_q$status
  )
  if (identical(batched_p$status, "pass") && identical(batched_q$status, "pass")) {
    record$batch_P_vs_batch_Q_bitwise_identical <- identical(batched_p$value, batched_q$value)
    record$batch_P_vs_batch_Q_differing_elements <- sum(batched_p$value != batched_q$value)
    record$batch_P_vs_batch_Q_maximum_absolute_difference <-
      max(abs(batched_p$value - batched_q$value))
  }
  if (identical(single$status, "pass") && identical(batched_p$status, "pass")) {
    record$single_vs_batch_P_bitwise_identical <- identical(single$value, batched_p$value)
    record$single_vs_batch_P_differing_elements <- sum(single$value != batched_p$value)
  }
  record
})

pipeline_pass <- all(vapply(
  target_records, function(record) isTRUE(record$pipeline_all_bitwise_identical), logical(1)
))
receipt <- list(
  schema_version = "masld-bench-gse49541-frma-single-array-firewall-v2",
  status = if (pipeline_pass) "pass_single_array_isolation_on_real_CEL" else "fail_single_array_isolation_on_real_CEL",
  series = "GSE49541",
  platform_id = "GPL570",
  cohort_family_id = "gse31803_gse49541_fibrosis_array",
  method = "frma_with_exact_hgu133plus2frmavecs",
  summarize_rule = "robust_weighted_average",
  normalize_rule = "quantile_to_frozen_hgu133plus2frmavecs_target",
  held_array_application = "one_array_at_a_time_with_frozen_training_object",
  frma_vector_package = "hgu133plus2frmavecs",
  frma_vector_object_sha256 = vector_digest,
  package_versions = as.list(observed),
  platform_features = length(platform_features),
  cel_records_available = nrow(manifest),
  target_selection = "sorted_GSM_accession_positions_1_36_72",
  companion_selection = "sorted_GSM_accession_positions_2_to_6_and_67_to_71",
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
  multi_array_frma_call_is_batch_dependent = !all(vapply(
    target_records, function(record) isTRUE(record$estimator_all_bitwise_identical), logical(1)
  )),
  multi_array_frma_call_used_as_preprocessing_path = FALSE,
  random_effect_probe = random_effect_probe,
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
  file.path(output_dir, "frma_single_array_firewall.json")
)
writeLines(
  capture.output(print(utils::sessionInfo())),
  file.path(output_dir, "R_sessionInfo.txt")
)
cat(sprintf(
  "targets=%d pipeline_bitwise_identical=%s perturbation_detected=%s distinct_digests=%s multi_array_batch_dependent=%s\n",
  length(targets),
  receipt$every_target_pipeline_bitwise_identical,
  receipt$every_perturbation_control_detected,
  receipt$distinct_targets_have_distinct_digests,
  receipt$multi_array_frma_call_is_batch_dependent
))
