#!/usr/bin/env Rscript

# Locate the stage of the frozen GPL570 fRMA path that depends on batch
# membership.
#
# model-check-085 proved on real CELs that the default fRMA configuration is not
# bit-identical between an array processed alone and the same array processed in
# a batch: about 70 percent of features move, by up to 9e-4 on the log2 scale.
# Batch order does not matter; batch membership does. This script finds where
# that dependence enters by holding one target array fixed and toggling the
# background, normalization, and summarization arguments one at a time.
#
# It also confirms the raw probe intensities of the target are identical between
# contexts, so any divergence is attributable to fRMA and not to file reading.
#
# Diagnostic only. No matrix is retained, no label is read, no model is fit.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: diagnose_frma_batch_dependence.R CEL_DIR MANIFEST_TSV AXIS_TSV OUTPUT_DIR")
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

axis <- utils::read.delim(
  axis_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
platform_features <- axis[["ID"]]

manifest <- utils::read.delim(
  manifest_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
manifest <- manifest[order(manifest$sample_accession), , drop = FALSE]
accessions <- manifest$sample_accession
cel_path <- function(accession) file.path(cel_dir, paste0(accession, ".CEL"))

target <- accessions[[1L]]
companion_p <- accessions[2L:6L]
companion_q <- accessions[67L:71L]

loader <- new.env(parent = emptyenv())
utils::data("hgu133plus2frmavecs", package = "hgu133plus2frmavecs", envir = loader)
frma_vectors <- loader[["hgu133plus2frmavecs"]]

## ------------- the raw intensities of the target must not depend on the batch
raw_of <- function(accession_set) {
  batch <- affy::ReadAffy(filenames = cel_path(accession_set))
  values <- affy::intensity(batch)
  as.numeric(values[, match(target, accession_set)])
}
raw_alone <- raw_of(target)
raw_in_p <- raw_of(c(target, companion_p))
raw_in_q <- raw_of(c(target, companion_q))
raw_identical <- identical(raw_alone, raw_in_p) && identical(raw_alone, raw_in_q)
rm(raw_in_p, raw_in_q)
invisible(gc(verbose = FALSE))

## ---------------------------------- one fRMA configuration across two contexts
summarize_target <- function(accession_set, background, normalize, summarize) {
  batch <- affy::ReadAffy(filenames = cel_path(accession_set))
  summary <- frma::frma(
    batch, background = background, normalize = normalize, summarize = summarize,
    input.vecs = frma_vectors, verbose = FALSE
  )
  values <- Biobase::exprs(summary)
  as.numeric(values[platform_features, match(target, accession_set)])
}

configurations <- list(
  list(id = "default_rma_quantile_rwa", background = "rma", normalize = "quantile", summarize = "robust_weighted_average"),
  list(id = "none_quantile_rwa",        background = "none", normalize = "quantile", summarize = "robust_weighted_average"),
  list(id = "rma_none_rwa",             background = "rma", normalize = "none",     summarize = "robust_weighted_average"),
  list(id = "none_none_rwa",            background = "none", normalize = "none",     summarize = "robust_weighted_average"),
  list(id = "rma_quantile_average",     background = "rma", normalize = "quantile", summarize = "average"),
  list(id = "rma_quantile_median",      background = "rma", normalize = "quantile", summarize = "median"),
  list(id = "rma_quantile_weighted_average", background = "rma", normalize = "quantile", summarize = "weighted_average")
)

records <- lapply(configurations, function(configuration) {
  attempt <- function(expression) {
    tryCatch(
      list(status = "pass", error = "not_applicable", value = expression),
      error = function(condition) list(status = "failed", error = conditionMessage(condition), value = NULL)
    )
  }
  alone <- attempt(summarize_target(target, configuration$background, configuration$normalize, configuration$summarize))
  if (!identical(alone$status, "pass")) {
    return(c(configuration, list(status = "failed", error = alone$error)))
  }
  in_p <- attempt(summarize_target(c(target, companion_p), configuration$background, configuration$normalize, configuration$summarize))
  in_q <- attempt(summarize_target(c(target, companion_q), configuration$background, configuration$normalize, configuration$summarize))
  if (!identical(in_p$status, "pass") || !identical(in_q$status, "pass")) {
    return(c(configuration, list(status = "failed", error = paste(in_p$error, in_q$error))))
  }
  reference <- alone$value
  c(configuration, list(
    status = "pass",
    error = "not_applicable",
    features = length(reference),
    finite = sum(is.finite(reference)),
    batch_P_bitwise_identical = identical(in_p$value, reference),
    batch_Q_bitwise_identical = identical(in_q$value, reference),
    batch_P_differing = sum(in_p$value != reference),
    batch_Q_differing = sum(in_q$value != reference),
    batch_P_max_absolute = max(abs(in_p$value - reference)),
    batch_Q_max_absolute = max(abs(in_q$value - reference)),
    batch_P_vs_Q_bitwise_identical = identical(in_p$value, in_q$value),
    reference_median = stats::median(reference),
    reference_iqr = stats::IQR(reference)
  ))
})

receipt <- list(
  schema_version = "masld-bench-frma-batch-dependence-diagnostic-v1",
  status = "pass_diagnostic_completed",
  series = "GSE49541",
  platform_id = "GPL570",
  target = target,
  companion_set_P = as.list(companion_p),
  companion_set_Q = as.list(companion_q),
  raw_intensities_identical_across_contexts = raw_identical,
  raw_probes = length(raw_alone),
  configurations = records,
  expression_values_exported = FALSE,
  labels_read = FALSE,
  model_training_activated = FALSE,
  sealed_outcomes_read = FALSE
)
writeLines(
  jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE, digits = NA),
  file.path(output_dir, "frma_batch_dependence_diagnostic.json")
)
# frmaAffyBatch is internal to the namespace; dump it only after the
# measurements above are already on disk.
tryCatch(
  writeLines(
    capture.output(print(utils::getFromNamespace("frmaAffyBatch", "frma"))),
    file.path(output_dir, "frmaAffyBatch_source.txt")
  ),
  error = function(condition) {
    writeLines(
      paste("frmaAffyBatch source unavailable:", conditionMessage(condition)),
      file.path(output_dir, "frmaAffyBatch_source.txt")
    )
  }
)
writeLines(
  capture.output(print(utils::sessionInfo())),
  file.path(output_dir, "R_sessionInfo.txt")
)
cat("diagnostic complete\n")
