#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
script_dir <- file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers")
Sys.setenv(HAC_ML_LIBRARY_ONLY = "1", HAC_ML_SCRIPT_DIR = script_dir)
source(file.path(script_dir, "21_tf_continuum.R"))

set.seed(20260817)
n <- 120L
axis <- rnorm(n)
fixture <- data.table(
  outcome_z = 0.7 * axis + rnorm(n, sd = 0.5),
  axis_raw = axis,
  fibrosis_stage = rep(0:3, length.out = n),
  inferred_sex = rep(c("F", "M"), length.out = n)
)
fit <- tf_fit_outcome_nonlinear(fixture, "fibrosis_stage")
stopifnot(fit$estimable, fit$beta > 0, fit$spline_df == 3L)

meta_fixture <- rbindlist(lapply(c("signature_pc1", "fixed_projection"), function(axis_id) {
  data.table(
    tf = rep(c("A", "B"), each = 2L),
    axis_id = axis_id,
    dataset = rep(c("GSE162694", "GSE213621"), 2L),
    estimable = TRUE,
    beta = c(0.5, 0.6, -0.3, -0.4),
    se = 0.1
  )
}))
meta <- tf_meta_models(meta_fixture, c("A", "B"),
                       c("signature_pc1", "fixed_projection"), 268L)
stopifnot(nrow(meta) == 4L, all(meta$direction_concordant),
          all(meta$bh_family_size == 268L))

source_text <- readLines(file.path(script_dir, "21_tf_continuum.R"), warn = FALSE)
stopifnot(
  any(grepl("1000L", source_text, fixed = TRUE)),
  any(grepl("tf_seed", source_text, fixed = TRUE)),
  any(grepl("signed_regulon_score", source_text, fixed = TRUE)),
  !any(grepl("install.packages|BiocManager::install|download.file", source_text))
)

set.seed(20260817)
mock_expression <- matrix(
  rnorm(12 * 20), nrow = 12,
  dimnames = list(paste0("G", 1:12), paste0("S", 1:20))
)
mock_metadata <- data.table(
  sample_id = colnames(mock_expression),
  dataset = rep(c("C1", "C2"), each = 10)
)
mock_network <- data.table(
  tf = rep(c("TF1", "TF2"), each = 6),
  target = rep(paste0("G", 1:6), 2),
  mor = rep(c(-1, 1), 6)
)
mock_scores <- tf_score_regulons(
  mock_expression, mock_metadata, mock_network, c("TF1", "TF2"),
  signature = character(), minimum_fraction = 0.8, minimum_targets = 5L
)
stopifnot(
  all(vapply(mock_scores$scores, function(x) any(is.finite(x)), logical(1))),
  all(vapply(mock_scores$scores, function(x) all(is.finite(x)), logical(1)))
)
cat("test_tf_contracts: PASS\n")
