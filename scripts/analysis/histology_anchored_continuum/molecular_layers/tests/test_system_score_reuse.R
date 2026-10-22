#!/usr/bin/env Rscript
# Competitive-null scoring must preserve whole-cohort normalization and each
# axis's eligible rows while reusing the outcome-independent collection medians.
suppressPackageStartupMessages(library(matrixStats))
root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
script_dir <- file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers")
local({
  source(file.path(script_dir, "lib_molecular_layers.R"), local = TRUE)
  expressions <- parse(file.path(script_dir, "73_stress_test_molecular_systems.R"))
  helper_names <- c("set_scores", "context_outcome", "fit_fast", "fixed_meta_z",
                    "statistic_for_set")
  for (expression in expressions) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% helper_names) {
      eval(expression)
    }
  }
  set.seed(20261021L)
  score_matrix <- matrix(rnorm(120L * 85L), nrow = 120L)
  collection_code <- rep(seq_len(7L), length.out = nrow(score_matrix))
  evaluation <- c("cohort_A", "cohort_B")
  axes_primary <- c("axis_A", "axis_B")
  contexts <- list()
  for (cohort_index in seq_along(evaluation)) {
    columns <- if (cohort_index == 1L) seq_len(40L) else 41:85
    stage <- rep(0:4, length.out = length(columns))
    sex <- rep(c("F", "M"), length.out = length(columns))
    for (axis_index in seq_along(axes_primary)) {
      # Different missing-axis and covariate rows must be removed only after
      # the whole-cohort score has been standardized.
      eligible <- setdiff(seq_along(columns), seq_len(2L + axis_index))
      nuisance <- model.matrix(~ factor(stage[eligible]) + factor(sex[eligible]))
      nuisance_qr <- qr(nuisance)
      axis <- ml_standardize(rnorm(length(eligible)))
      residual <- as.numeric(qr.resid(nuisance_qr, axis))
      residual <- residual - mean(residual)
      contexts[[paste(axes_primary[[axis_index]], evaluation[[cohort_index]], sep = "|")]] <- list(
        donor_columns = columns, analysis_positions = eligible,
        nuisance_qr = nuisance_qr, x_residual = residual,
        sxx = sum(residual^2), residual_df = length(eligible) - ncol(nuisance) - 1L
      )
    }
  }
  reference <- function(node_indices) {
    vapply(axes_primary, function(axis_name) {
      fits <- lapply(evaluation, function(cohort) {
        context <- contexts[[paste(axis_name, cohort, sep = "|")]]
        fit_fast(context, set_scores(node_indices, context$donor_columns))
      })
      abs(fixed_meta_z(vapply(fits, `[[`, numeric(1), "beta"),
                       vapply(fits, `[[`, numeric(1), "se")))
    }, numeric(1))
  }
  for (size in c(1L, 7L, 23L, 120L)) {
    for (draw in seq_len(5L)) {
      nodes <- sample.int(nrow(score_matrix), size)
      saved_contexts <- serialize(contexts, NULL)
      before_rng <- .Random.seed
      expected <- reference(nodes)
      stopifnot(identical(.Random.seed, before_rng))
      actual <- statistic_for_set(nodes)
      stopifnot(identical(actual, expected), identical(.Random.seed, before_rng),
                identical(serialize(contexts, NULL), saved_contexts))
    }
  }
  # Verify that the expensive medians are evaluated once per cohort.
  scoring_calls <- 0L
  score_original <- set_scores
  set_scores <- function(...) {
    scoring_calls <<- scoring_calls + 1L
    score_original(...)
  }
  invisible(statistic_for_set(seq_len(23L)))
  stopifnot(scoring_calls == length(evaluation))
})
cat("SYSTEM_SCORE_REUSE_TESTS_PASS\n")
