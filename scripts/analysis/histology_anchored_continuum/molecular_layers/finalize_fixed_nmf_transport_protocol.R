#!/usr/bin/env Rscript
# Bind the prepared protocol to its separate normalization checks and unavailable axes.
suppressPackageStartupMessages({library(data.table); library(NMF); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: protocol_review_directory transport_check_directory")
path <- file.path(args[1L], "revised_fixed_nmf_projection_protocol_with_transport_checks_r2.rds")
if (file.exists(path)) stop("Refusing to overwrite a reviewed protocol")
protocol <- readRDS(file.path(args[1L], "revised_fixed_nmf_projection_protocol.rds"))
preparation <- fromJSON(file.path(args[1L], "projection_preparation.json"))
checks <- fread(file.path(args[2L], "native_normalization_transport_checks.tsv"))
native_review <- fromJSON(file.path(args[2L], "native_transport_review.json"))
stopifnot(nrow(checks) == 20L, !anyDuplicated(checks[, paste(k, program, variant)]),
  all(checks$all_native_scores_finite_nonnegative), !native_review$corrected_scores_produced)
expected <- with(checks, historical_global_spearman >= .99 & historical_min_cohort_spearman >= .9 &
  retained_basis_mass >= .8 & retained_basis_rank == k)
stopifnot(identical(expected, checks$transport_reconstruction_eligible))
original <- readRDS(protocol$historical_preprocessing$source_paths$cache)
for (k in c("4", "6")) {
  full <- NMF::basis(original$nmf_results[[k]])
  stopifnot(identical(protocol$basis[[k]], full[protocol$basis_feature_mask$retained, , drop = FALSE]))
}
protocol$historical_normalization_transport_checks <- checks
protocol$annotation_correspondence_sha256 <- preparation$source_sha256$annotation
unavailable_axes <- function(variant_name) {
  selected <- checks[variant == variant_name & !transport_reconstruction_eligible]
  if (!nrow(selected)) return(character())
  paste0("k", selected$k, "-", selected$program)
}
protocol$stage_unavailable_axes <- unavailable_axes("full_stage")
protocol$continuum_unavailable_axes <- unavailable_axes("signature_excluded_continuum")
stopifnot(length(protocol$stage_unavailable_axes) == sum(!checks[variant == "full_stage", transport_reconstruction_eligible]),
          length(protocol$continuum_unavailable_axes) == sum(!checks[variant == "signature_excluded_continuum", transport_reconstruction_eligible]))
protocol$unavailable_score_policy <- "Keep every axis and unit row; unavailable axes have NA scores and explicit failed-reconstruction reasons, never zero or outcome-derived replacement. Do not submit partially unavailable scores to the historical finite-score stage input contract."
protocol$primary_inference_ready <- FALSE
protocol$corrected_count_projection_ready <- FALSE
protocol$transport_check_path <- normalizePath(args[2L])
saveRDS(protocol, path)
report <- list(status = "prepared_revised_sensitivity_protocol_with_unavailable_axes",
  annotation_unresolved = preparation$annotation_unresolved, original_weights_unchanged = TRUE,
  corrected_scores_produced = FALSE, original_axes_replaced = FALSE,
  stage_unavailable_axes = protocol$stage_unavailable_axes,
  continuum_unavailable_axes = protocol$continuum_unavailable_axes,
  protocol_path = normalizePath(path), fixed_criteria_relaxed = FALSE)
writeLines(toJSON(report, auto_unbox = TRUE, pretty = TRUE), file.path(args[1L], "final_protocol_review_r2.json"))
cat(toJSON(report, auto_unbox = TRUE), "\n")
