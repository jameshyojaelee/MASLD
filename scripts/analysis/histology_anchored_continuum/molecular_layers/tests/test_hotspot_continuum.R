#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
script_dir <- file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers")
Sys.setenv(
  HAC_ML_TEST_MODE = "1", HAC_ML_PERMUTATIONS = "100",
  HAC_ML_SCRIPT_DIR = script_dir
)
source(file.path(script_dir, "31_hotspot_continuum.R"))

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

set.seed(20260817)
n <- 160L
stage <- rep(0:3, each = n / 4)
sex <- rep(c("F", "M"), length.out = n)
axis <- as.numeric(scale(stage + rnorm(n)))
fixture <- data.table(
  outcome_z = 0.7 * axis + 0.2 * stage + rnorm(n, sd = 0.7),
  axis_raw = axis,
  fibrosis_stage = stage,
  inferred_sex = sex
)
fit <- hotspot_fit_continuum(fixture, "fibrosis_stage")
assert(fit$estimable && fit$beta > 0,
       "Hotspot stage/sex model lost a positive within-stage continuum effect")

permutation_a <- hotspot_within_stage_permutation(fixture, 100L, 8123L)
permutation_b <- hotspot_within_stage_permutation(fixture, 100L, 8123L)
assert(identical(permutation_a, permutation_b),
       "Within-stage Hotspot permutations are not deterministic")
assert(permutation_a$permutation_replicates == 100L &&
         permutation_a$n_finite_permutations == 100L,
       "Permutation accounting is incomplete")

adjusted <- ml_complete_bh(c(0.001, rep(NA_real_, 116L)), 117L)
assert(abs(adjusted[[1L]] - 0.117) < 1e-12,
       "Hotspot BH adjustment did not retain the complete n=117 family")

source_text <- readLines(file.path(script_dir, "31_hotspot_continuum.R"), warn = FALSE)
assert(any(grepl("fibrosis_stage_by_inferred_sex", source_text, fixed = TRUE)),
       "Permutation strata contract is absent")
assert(any(grepl("signature_excluded", source_text, fixed = TRUE)),
       "Signature-gene exclusion contract is absent")
assert(any(grepl("10000L", source_text, fixed = TRUE)),
       "Production 10,000-permutation gate is absent")
cat("test_hotspot_continuum: PASS\n")
