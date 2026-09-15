#!/usr/bin/env Rscript

# Verify that the repaired matrix writer round-trips IEEE-754 doubles exactly,
# using the real GSE49541 values preserved in the failed stage of job 21109357.
#
# utils::write.table converts doubles at 15 significant digits. A double needs
# 17 to survive a text round trip, so the first write of the platform matrix was
# a lossy representation of what fRMA actually produced. The bitwise assertion
# in the summarizer caught it; this script confirms the repair on real values at
# real magnitudes rather than on synthetic ones, and confirms that the old
# behaviour genuinely fails, so the test is not vacuous.
#
# Diagnostic only. Reads no label, fits no model, publishes no matrix.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: verify_matrix_text_roundtrip.R LOSSY_MATRIX_TSV OUTPUT_DIR")
source_matrix <- args[[1L]]
output_dir <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(jsonlite))

# The preserved matrix cannot be used as test data directly. Its values were
# already written at 15 significant digits, so each one is the nearest double to
# a 15-digit decimal. Re-writing those at 15 digits reproduces the same decimal
# and therefore the same double, and the test passes no matter what the writer
# does. The first run of this script made exactly that mistake and reported all
# three writers as lossless.
#
# The values are therefore perturbed by a sub-ULP-scale multiplier to restore
# full mantissa entropy, giving generic doubles at genuine fRMA magnitudes -
# what the summarizer actually held in memory. The default writer is then
# REQUIRED to fail, which is what makes a pass by the repaired writer meaningful.
frame <- utils::read.delim(
  source_matrix, sep = "\t", quote = "", check.names = FALSE, row.names = 1L
)
values <- as.matrix(frame)
rm(frame)
invisible(gc(verbose = FALSE))

set.seed(20260825L)
values <- values * (1 + stats::runif(length(values), 1e-13, 1e-12))
fifteen_digit_clean <- sum(
  as.numeric(sprintf("%.15g", values)) == as.numeric(values)
)

round_trip <- function(values, digits, path, use_sprintf) {
  if (use_sprintf) {
    text <- matrix(
      sprintf(paste0("%.", digits, "g"), values),
      nrow = nrow(values), ncol = ncol(values)
    )
    out <- data.frame(
      id = rownames(values), text, check.names = FALSE, stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      id = rownames(values), values, check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  names(out) <- c("platform_feature_id", colnames(values))
  utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  restored <- as.matrix(utils::read.delim(
    path, sep = "\t", quote = "", check.names = FALSE, row.names = 1L
  ))
  identical_bits <- identical(dim(restored), dim(values)) &&
    identical(rownames(restored), rownames(values)) &&
    identical(colnames(restored), colnames(values)) &&
    identical(as.numeric(restored), as.numeric(values))
  differing <- sum(as.numeric(restored) != as.numeric(values))
  maximum <- max(abs(as.numeric(restored) - as.numeric(values)))
  unlink(path)
  list(
    bitwise_identical = identical_bits,
    differing_elements = differing,
    maximum_absolute_difference = maximum,
    bytes_written = NA_real_
  )
}

default_writer <- round_trip(values, NA, file.path(output_dir, "default.tsv"), FALSE)
repaired_writer <- round_trip(values, 17L, file.path(output_dir, "repaired.tsv"), TRUE)
fifteen_digits <- round_trip(values, 15L, file.path(output_dir, "fifteen.tsv"), TRUE)

receipt <- list(
  schema_version = "masld-bench-matrix-text-roundtrip-v1",
  status = if (isTRUE(repaired_writer$bitwise_identical) &&
                !isTRUE(default_writer$bitwise_identical) &&
                !isTRUE(fifteen_digits$bitwise_identical))
    "pass_repaired_writer_round_trips_and_default_does_not"
  else "fail_round_trip_evidence_is_inconclusive",
  source_matrix = source_matrix,
  test_values_perturbed_to_restore_mantissa_entropy = TRUE,
  values_already_clean_at_15_significant_digits = fifteen_digit_clean,
  negative_control_required = "default and 15-digit writers must both fail",
  features = nrow(values),
  arrays = ncol(values),
  values_tested = length(values),
  default_write_table = default_writer,
  fifteen_significant_digits = fifteen_digits,
  seventeen_significant_digits = repaired_writer,
  finding = paste(
    "utils::write.table serialises doubles at 15 significant digits.",
    "IEEE-754 doubles need 17. A bitwise round-trip assertion after every",
    "numeric matrix write is the cheap way to catch this."
  ),
  labels_read = FALSE,
  model_training_activated = FALSE
)
writeLines(
  jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE, digits = NA),
  file.path(output_dir, "matrix_text_roundtrip.json")
)
cat(sprintf(
  "default_bitwise=%s 15g_bitwise=%s 17g_bitwise=%s values=%d\n",
  default_writer$bitwise_identical, fifteen_digits$bitwise_identical,
  repaired_writer$bitwise_identical, length(values)
))
