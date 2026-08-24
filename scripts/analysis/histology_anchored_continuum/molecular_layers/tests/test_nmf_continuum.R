#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
script_dir <- file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers")
Sys.setenv(HAC_ML_TEST_MODE = "1", HAC_ML_SCRIPT_DIR = script_dir)
source(file.path(script_dir, "30_nmf_continuum.R"))

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

set.seed(20260817)
n <- 120L
axis <- rnorm(n)
fixture <- data.table(
  outcome_z = 0.8 * axis + rnorm(n, sd = 0.5),
  axis_raw = axis,
  fibrosis_stage = rep(0:3, length.out = n),
  inferred_sex = rep(c("F", "M"), length.out = n)
)
fit <- nmf_fit_continuum(fixture, "fibrosis_stage")
assert(fit$estimable && fit$beta > 0, "NMF stage/sex model lost a positive continuum effect")
shape <- nmf_adjusted_shape(fixture, "fibrosis_stage")
assert(nrow(shape) == 101L && identical(shape$percentile, seq(0, 1, length.out = 101L)),
       "NMF spline prediction grid is not the frozen 101 percentiles")

if (requireNamespace("nnls", quietly = TRUE)) {
  basis <- matrix(runif(240L, 0.1, 2), nrow = 60L, ncol = 4L)
  colnames(basis) <- paste0("P", 1:4)
  coefficients <- matrix(runif(32L, 0, 1), nrow = 4L)
  expression <- basis %*% coefficients
  reconstructed <- nmf_nnls_coefficients(basis, expression)
  assert(max(abs(reconstructed - coefficients)) < 1e-7,
         "Fixed-W NNLS does not recover exact synthetic coefficients")
}

adjusted <- ml_complete_bh(c(0.001, rep(NA_real_, 9L)), 10L)
assert(abs(adjusted[[1L]] - 0.01) < 1e-12,
       "NMF BH adjustment did not retain the complete n=10 family")

source_text <- readLines(file.path(script_dir, "30_nmf_continuum.R"), warn = FALSE)
assert(!any(grepl("A3_axis|HAC_NMF_AXIS_PATH", source_text, fixed = FALSE)),
       "The forbidden target-oriented A3 composite entered NMF inference")
assert(any(grepl("reconstruction_gate_pass", source_text, fixed = TRUE)),
       "Fixed-W reconstruction gates are absent")
cat("test_nmf_continuum: PASS\n")
