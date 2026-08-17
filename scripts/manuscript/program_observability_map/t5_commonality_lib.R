# T5 composition coupling: variance-decomposition primitives.
#
# WHAT THIS COMPUTES AND WHAT IT DOES NOT
# For one program score y, one cohort, and three blocks of predictors
#   C = nuisance covariates always in the model (intercept, sex, NAS)
#   A = the fibrosis term
#   B = the accepted composition block (3 within-cohort CLR principal components)
# this returns the commonality partition of the model R-squared:
#
#   R2(C+A+B) = R2(C) + unique_A + unique_B + common(A,B)
#
# with
#   unique_A = R2(C+A+B) - R2(C+B)
#   unique_B = R2(C+A+B) - R2(C+A)
#   common   = R2(C+A) + R2(C+B) - R2(C+A+B) - R2(C)
#
# The identity is exact by construction and t5_tests.R checks it to 1e-12.
#
# `common` is signed. A negative value is suppression: composition and fibrosis
# each explain more of the program in the other's presence than alone. It is a
# real and interpretable configuration, not an error, and it is NEVER clamped to
# zero here. Clamping would silently convert suppression into "no coupling".
#
# THE CEILING. The composition proportions are deconvolved from the same logCPM
# matrix the program scores are computed from. A program built out of a lineage's
# marker genes therefore shares variance with that lineage's estimated abundance
# partly by construction. `common` is consequently an UPPER BOUND on how much of
# a program's fibrosis association could be composition, and nothing in this file
# can separate the mechanical part of that bound from the biological part. Report
# coupling, never attribution.

suppressPackageStartupMessages({
  library(data.table)
})

# Residual sums of squares for every column of Y against one shared design.
# Returns NULL when the design is rank deficient, so a caller can refuse rather
# than silently fit a projected version of the model it asked for.
t5_rss <- function(Y, X) {
  qrX <- qr(X)
  if (qrX$rank < ncol(X) || nrow(X) <= ncol(X)) return(NULL)
  E <- qr.resid(qrX, Y)
  colSums(E * E)
}

# The parts of the decomposition that do not involve the composition block. They
# are invariant under a composition permutation, so precomputing them makes a
# 3,000-draw null affordable without changing any number it produces.
t5_baseline_parts <- function(Y, C, A) {
  centred <- sweep(Y, 2L, colMeans(Y), "-")
  tss <- colSums(centred * centred)
  if (any(!is.finite(tss)) || any(tss <= 0)) {
    fail("A program score is constant within a cohort; its R-squared is undefined")
  }
  list(tss = tss, rss_C = t5_rss(Y, C), rss_CA = t5_rss(Y, cbind(C, A)))
}

# One cohort, all programs at once. Y is samples x programs.
t5_commonality_cohort <- function(Y, C, A, B, parts = NULL) {
  stopifnot(is.matrix(Y), nrow(Y) == nrow(C), nrow(Y) == nrow(A), nrow(Y) == nrow(B))
  n <- nrow(Y)
  if (is.null(parts)) parts <- t5_baseline_parts(Y, C, A)
  tss <- parts$tss
  rss_C <- parts$rss_C
  rss_CA <- parts$rss_CA
  rss_CB <- t5_rss(Y, cbind(C, B))
  rss_CAB <- t5_rss(Y, cbind(C, A, B))
  if (is.null(rss_C) || is.null(rss_CA) || is.null(rss_CB) || is.null(rss_CAB)) {
    return(NULL)
  }
  r2 <- function(rss) 1 - rss / tss
  R2_C <- r2(rss_C)
  R2_CA <- r2(rss_CA)
  R2_CB <- r2(rss_CB)
  R2_CAB <- r2(rss_CAB)

  df_A <- ncol(A)
  df_B <- ncol(B)
  df_error <- n - (ncol(C) + df_A + df_B)
  mse <- rss_CAB / df_error
  f_unique_A <- ((rss_CB - rss_CAB) / df_A) / mse
  f_unique_B <- ((rss_CA - rss_CAB) / df_B) / mse

  list(
    n = n,
    r2_covariates = R2_C,
    r2_covariates_fibrosis = R2_CA,
    r2_covariates_composition = R2_CB,
    r2_full = R2_CAB,
    unique_fibrosis = R2_CAB - R2_CB,
    unique_composition = R2_CAB - R2_CA,
    common = R2_CA + R2_CB - R2_CAB - R2_C,
    fibrosis_incremental = R2_CA - R2_C,
    f_unique_fibrosis = f_unique_A,
    f_unique_composition = f_unique_B,
    df_unique_fibrosis = df_A,
    df_unique_composition = df_B,
    df_error = df_error,
    p_unique_fibrosis = stats::pf(f_unique_A, df_A, df_error, lower.tail = FALSE),
    p_unique_composition = stats::pf(f_unique_B, df_B, df_error, lower.tail = FALSE)
  )
}

# Donor-weighted aggregation of per-cohort components. R-squared is a
# within-cohort quantity (scores are standardized within cohort and fibrosis is
# centered within cohort), so the cohorts are combined by weighting each one by
# its donor count. Nothing is pooled across cohorts at the sample level.
t5_aggregate_cohorts <- function(per_cohort) {
  weights <- vapply(per_cohort, function(x) x$n, numeric(1))
  weights <- weights / sum(weights)
  combine <- function(field) {
    Reduce(`+`, Map(function(x, w) x[[field]] * w, per_cohort, weights))
  }
  fields <- c("r2_covariates", "r2_covariates_fibrosis", "r2_covariates_composition",
              "r2_full", "unique_fibrosis", "unique_composition", "common",
              "fibrosis_incremental")
  out <- lapply(fields, combine)
  names(out) <- fields
  out$n_cohorts <- length(per_cohort)
  out$n_donors <- sum(vapply(per_cohort, function(x) x$n, numeric(1)))
  out
}

# Cohorts are disjoint donor sets, so their partial-F tests are independent and
# Fisher's method is exact under the global null. This tests "no unique
# contribution in ANY cohort", which is a different and weaker estimand than the
# meta-analysed common effect the axis map reports; both are written out.
t5_fisher_combine <- function(p_matrix) {
  stat <- -2 * rowSums(log(pmax(p_matrix, .Machine$double.xmin)))
  df <- 2 * ncol(p_matrix)
  list(statistic = stat, df = df,
       p_value = stats::pchisq(stat, df = df, lower.tail = FALSE))
}

# The reported coupling quantity.
#
#   raw fraction       = common / (R2(C+A) - R2(C))
#   corrected fraction = (common - null mean of common) / (R2(C+A) - R2(C))
#
# The denominator is the fibrosis-associated increment over the nuisance
# covariates. It contains no composition term, so it is identical in the observed
# data and under a composition permutation, which is what makes the correction
# well defined.
#
# The correction exists because unadjusted R-squared rises when parameters are
# added even if they carry no information: three permuted composition columns
# manufacture a positive `common` on their own. The corrected fraction subtracts
# that floor and is the honest headline. The raw fraction is kept because only
# the raw terms satisfy the exact additivity identity.
t5_coupling_fraction <- function(common, fibrosis_incremental, common_null_mean = 0) {
  usable <- is.finite(common) & is.finite(fibrosis_incremental) &
    fibrosis_incremental > 0
  ifelse(usable, (common - common_null_mean) / fibrosis_incremental, NA_real_)
}

# Joint within-cohort permutation of a covariate block. All columns move with the
# same permutation so their mutual correlation is preserved and only the link to
# the donor's own expression is broken.
t5_permute_block_within <- function(block, group, strata = NULL) {
  out <- block
  key <- if (is.null(strata)) as.character(group) else paste(group, strata, sep = "\r")
  for (level in unique(key)) {
    rows <- which(key == level)
    if (length(rows) < 2L) next
    out[rows, ] <- block[sample(rows), , drop = FALSE]
  }
  out
}
