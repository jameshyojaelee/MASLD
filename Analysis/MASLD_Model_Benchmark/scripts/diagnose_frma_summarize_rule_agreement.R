#!/usr/bin/env Rscript

# Measure how far fRMA's batch-invariant summarization rules sit from its
# default, on arrays processed alone.
#
# model-check-087 showed that only summarize = "robust_weighted_average" carries
# batch dependence, and that "weighted_average", "average", and "median" are
# bit-identical whatever the batch. The remaining question is whether pinning a
# batch-invariant rule is a specification of the frozen method or a change to it,
# so this compares the rules on the same single arrays.
#
# Diagnostic only. No matrix is retained, no label is read, no model is fit.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: diagnose_frma_summarize_rule_agreement.R CEL_DIR MANIFEST_TSV AXIS_TSV OUTPUT_DIR")
}
cel_dir <- args[[1L]]; manifest_path <- args[[2L]]
axis_path <- args[[3L]]; output_dir <- args[[4L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages({ library(affy); library(frma); library(jsonlite) })

axis <- utils::read.delim(axis_path, sep = "\t", quote = "", check.names = FALSE,
                          colClasses = "character", na.strings = character(0))
platform_features <- axis[["ID"]]
manifest <- utils::read.delim(manifest_path, sep = "\t", quote = "", check.names = FALSE,
                              colClasses = "character", na.strings = character(0))
manifest <- manifest[order(manifest$sample_accession), , drop = FALSE]
accessions <- manifest$sample_accession[c(1L, 36L, 72L)]

loader <- new.env(parent = emptyenv())
utils::data("hgu133plus2frmavecs", package = "hgu133plus2frmavecs", envir = loader)
frma_vectors <- loader[["hgu133plus2frmavecs"]]

summarize_alone <- function(accession, rule) {
  batch <- affy::ReadAffy(filenames = file.path(cel_dir, paste0(accession, ".CEL")))
  values <- Biobase::exprs(frma::frma(batch, summarize = rule,
                                      input.vecs = frma_vectors, verbose = FALSE))
  as.numeric(values[platform_features, 1L])
}

rules <- c("weighted_average", "average", "median")
records <- lapply(accessions, function(accession) {
  reference <- summarize_alone(accession, "robust_weighted_average")
  comparisons <- lapply(rules, function(rule) {
    alternative <- summarize_alone(accession, rule)
    difference <- alternative - reference
    list(
      rule = rule,
      pearson_r = stats::cor(alternative, reference),
      spearman_rho = stats::cor(alternative, reference, method = "spearman"),
      median_absolute_difference = stats::median(abs(difference)),
      q95_absolute_difference = unname(stats::quantile(abs(difference), 0.95)),
      max_absolute_difference = max(abs(difference)),
      features_differing_by_more_than_0_01 = sum(abs(difference) > 0.01),
      features_differing_by_more_than_0_10 = sum(abs(difference) > 0.10)
    )
  })
  list(
    sample_accession = accession,
    reference_rule = "robust_weighted_average",
    reference_median = stats::median(reference),
    reference_iqr = stats::IQR(reference),
    comparisons = comparisons
  )
})

writeLines(
  jsonlite::toJSON(list(
    schema_version = "masld-bench-frma-summarize-rule-agreement-v1",
    status = "pass_diagnostic_completed",
    series = "GSE49541", platform_id = "GPL570",
    arrays = as.list(accessions), rules_compared = as.list(rules),
    per_array = records,
    expression_values_exported = FALSE, labels_read = FALSE,
    model_training_activated = FALSE, sealed_outcomes_read = FALSE
  ), auto_unbox = TRUE, pretty = TRUE, digits = NA),
  file.path(output_dir, "frma_summarize_rule_agreement.json")
)
cat("agreement diagnostic complete\n")
