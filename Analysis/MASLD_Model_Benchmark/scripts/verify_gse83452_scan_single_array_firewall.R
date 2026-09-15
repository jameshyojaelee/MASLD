#!/usr/bin/env Rscript

# Put SCAN through the same four pipeline contexts as oligo, so the question
# "does SCAN pool across arrays?" is answered with evidence rather than left
# either assumed safe or assumed unsafe.
#
# What is already known: SCAN runs on a single array once foreach has a
# sequential backend registered (its %dopar% failure was a SCAN.UPC packaging
# defect, not a missing dependency), and a repeat call is bit-identical. That is
# determinism, NOT a demonstration that it avoids cross-array pooling.
#
# The batching route cannot answer it either: a multi-array SCAN call errors
# inside SCAN.UPC. So the only route left is the pipeline route, which is also
# the one that matters: summarize a target alone, and again with different work
# before and after it in the same session, and require bit-identity.
#
# A batched call is still attempted and its outcome recorded, because "SCAN
# cannot be batched at all" is itself a safety-relevant property worth freezing.
#
# 1e-12 perturbation control and blind positional target selection, unchanged.
# Reads no label, no GEO series matrix, exports no expression value, fits no model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: verify_gse83452_scan_single_array_firewall.R CEL_DIR MANIFEST_TSV OUTPUT_DIR")
}
cel_dir <- args[[1L]]; manifest_path <- args[[2L]]; output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages({
  library(SCAN.UPC); library(pd.hugene.2.0.st); library(foreach); library(jsonlite)
})
foreach::registerDoSEQ()

EXPECTED <- c("SCAN.UPC" = "2.52.0", "pd.hugene.2.0.st" = "3.14.1", "oligo" = "1.74.0")
observed <- vapply(names(EXPECTED), function(p) as.character(packageVersion(p)), character(1))
if (!identical(unname(observed), unname(EXPECTED))) {
  stop(sprintf("SCAN package versions differ: %s",
               paste(sprintf("%s=%s", names(observed), observed), collapse = ",")))
}

manifest <- utils::read.delim(manifest_path, sep = "\t", quote = "", check.names = FALSE,
                              colClasses = "character", na.strings = character(0))
manifest <- manifest[manifest$series == "GSE83452", , drop = FALSE]
manifest <- manifest[order(manifest$sample_accession), , drop = FALSE]
accessions <- manifest$sample_accession
cel_path <- function(a) file.path(cel_dir, paste0(a, ".CEL"))

targets <- accessions[c(1L, 116L, 231L)]
companion_p <- accessions[2L:4L]
companion_q <- accessions[228L:230L]

feature_axis <- NULL
digest_doubles <- function(values) {
  path <- tempfile(fileext = ".f8"); on.exit(unlink(path), add = TRUE)
  con <- file(path, "wb"); writeBin(as.double(values), con, size = 8L, endian = "little"); close(con)
  substr(system2("sha256sum", shQuote(path), stdout = TRUE), 1L, 64L)
}
scan_one <- function(accession) {
  values <- Biobase::exprs(SCAN.UPC::SCAN(cel_path(accession), outFilePath = NA, verbose = FALSE))
  if (ncol(values) != 1L) stop(sprintf("%s was not summarized as a single array", accession))
  if (is.null(feature_axis)) feature_axis <<- rownames(values)
  else if (!identical(rownames(values), feature_axis)) stop("SCAN feature axis drifted")
  out <- as.numeric(values[feature_axis, 1L]); rm(values); invisible(gc(verbose = FALSE)); out
}

pipeline_contexts <- list(
  solo_only = function(t) scan_one(t),
  solo_after_P = function(t) { for (c in companion_p) invisible(scan_one(c)); scan_one(t) },
  solo_before_P = function(t) { v <- scan_one(t); for (c in companion_p) invisible(scan_one(c)); v },
  solo_interleaved_Q = function(t) {
    invisible(scan_one(companion_q[[1L]])); v <- scan_one(t)
    for (c in companion_q[2L:length(companion_q)]) invisible(scan_one(c)); v
  }
)

records <- list(); alone <- character(0)
for (target in targets) {
  vals <- lapply(pipeline_contexts, function(f) f(target))
  ref <- vals[["solo_only"]]
  if (!all(is.finite(ref))) stop("SCAN summary carries a nonfinite value")
  digests <- vapply(vals, digest_doubles, character(1))
  perturbed <- ref; perturbed[[1L]] <- perturbed[[1L]] + 1e-12
  records[[target]] <- list(
    sample_accession = target,
    pipeline_digests = as.list(digests),
    pipeline_all_digests_identical = length(unique(unname(digests))) == 1L,
    pipeline_all_bitwise_identical = all(vapply(vals, function(v) identical(v, ref), logical(1))),
    pipeline_differing_elements = lapply(vals, function(v) sum(v != ref)),
    pipeline_maximum_absolute_difference = lapply(vals, function(v) max(abs(v - ref))),
    pipeline_arrays_per_summarization_call = 1L,
    perturbation_control_digest_changed = !identical(digest_doubles(perturbed), digests[["solo_only"]]),
    perturbation_control_magnitude = 1e-12,
    features = length(ref), finite_values = sum(is.finite(ref))
  )
  alone <- c(alone, unname(digests[["solo_only"]]))
  rm(vals); invisible(gc(verbose = FALSE))
}

batched <- tryCatch({
  v <- Biobase::exprs(SCAN.UPC::SCAN(cel_path(c(targets[[1L]], companion_p)),
                                     outFilePath = NA, verbose = FALSE))
  list(status = "pass", columns = ncol(v))
}, error = function(e) list(status = "failed", error = conditionMessage(e)))

pass <- all(vapply(records, function(r) isTRUE(r$pipeline_all_bitwise_identical), logical(1)))
receipt <- list(
  schema_version = "masld-bench-gse83452-scan-single-array-firewall-v1",
  status = if (pass) "pass_single_array_isolation_on_real_CEL" else "fail_single_array_isolation_on_real_CEL",
  series = "GSE83452", platform_id = "GPL16686",
  method = "SCAN_with_pd.hugene.2.0.st",
  held_array_application = "one_array_at_a_time_with_frozen_training_object",
  foreach_sequential_backend_registered = TRUE,
  package_versions = as.list(observed),
  scan_features = if (is.null(feature_axis)) NA_integer_ else length(feature_axis),
  cel_records_available = nrow(manifest),
  target_selection = "sorted_GSM_accession_positions_1_116_231",
  selection_used_labels = FALSE, selection_used_intensity = FALSE,
  targets = as.list(targets),
  per_target = unname(records),
  every_target_pipeline_bitwise_identical = pass,
  every_perturbation_control_detected = all(vapply(
    records, function(r) isTRUE(r$perturbation_control_digest_changed), logical(1))),
  distinct_targets_have_distinct_digests = length(unique(alone)) == length(alone),
  multi_array_call = batched,
  multi_array_call_is_possible = identical(batched$status, "pass"),
  scan_pools_across_arrays_when_batched = if (identical(batched$status, "pass"))
    "measured_separately_not_in_this_run" else "not_applicable_multi_array_call_errors",
  gpl570_evidence_carried_over = FALSE,
  all_sample_RMA_run = FALSE, across_array_quantile_normalization_run = FALSE,
  GEO_series_matrix_read = FALSE, expression_values_exported = FALSE,
  labels_read = FALSE, model_training_activated = FALSE, sealed_outcomes_read = FALSE
)
writeLines(jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE, digits = NA),
           file.path(output_dir, "scan_single_array_firewall.json"))
writeLines(capture.output(print(utils::sessionInfo())), file.path(output_dir, "R_sessionInfo.txt"))
cat(sprintf("scan_features=%s pipeline_bitwise=%s perturbation=%s distinct=%s batched=%s\n",
            length(feature_axis), receipt$every_target_pipeline_bitwise_identical,
            receipt$every_perturbation_control_detected,
            receipt$distinct_targets_have_distinct_digests, batched$status))
