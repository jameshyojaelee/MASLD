#!/usr/bin/env Rscript

# T5 unit tests. Synthetic data with a decomposition known in closed form, run
# before the real substrate is touched. The run script refuses to start unless
# this file exits 0.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(script_dir, "t5_commonality_lib.R"))

passed <- 0L
check <- function(condition, label) {
  if (!isTRUE(condition)) fail("T5 TEST FAILED: ", label)
  passed <<- passed + 1L
  message("  ok  ", label)
}

set.seed(20260812L)

# --------------------------------------------------------------- 1. additivity
# The identity must hold on arbitrary correlated data, not only on a tidy case.
n <- 300L
C <- cbind(intercept = 1, sex = rbinom(n, 1L, 0.4), nas = rnorm(n))
A <- matrix(0.5 * C[, "nas"] + rnorm(n), ncol = 1L, dimnames = list(NULL, "fib"))
B <- matrix(rnorm(n * 3L), ncol = 3L, dimnames = list(NULL, paste0("pc", 1:3)))
B[, 1L] <- B[, 1L] + 0.6 * A[, 1L]
Y <- cbind(
  y1 = 0.8 * A[, 1L] + 0.4 * B[, 1L] + rnorm(n),
  y2 = 0.2 * C[, "nas"] + rnorm(n),
  y3 = rnorm(n)
)
fit <- t5_commonality_cohort(Y, C, A, B)
residual <- with(fit, r2_covariates + unique_fibrosis + unique_composition +
                   common - r2_full)
check(max(abs(residual)) < 1e-12,
      sprintf("components sum to the full R-squared (max |error| %.3g)", max(abs(residual))))

# --------------------------------------------------- 2. closed-form known case
# Orthonormal A and B, both orthogonal to the intercept and to each other, with
# an error term orthogonal to all of them. Every component is known exactly.
orth <- qr.Q(qr(cbind(1, matrix(rnorm(n * 3L), ncol = 3L))))
C2 <- orth[, 1L, drop = FALSE]
A2 <- orth[, 2L, drop = FALSE]
B2 <- orth[, 3:4, drop = FALSE]
e <- qr.resid(qr(orth), rnorm(n))
a <- 3; b <- 2
y <- a * A2[, 1L] + b * B2[, 1L] + e
Y2 <- matrix(y, ncol = 1L, dimnames = list(NULL, "y"))
fit2 <- t5_commonality_cohort(Y2, C2, A2, B2)
tss <- sum((y - mean(y))^2)
check(abs(fit2$unique_fibrosis - a^2 / tss) < 1e-12, "unique fibrosis matches closed form")
check(abs(fit2$unique_composition - b^2 / tss) < 1e-12, "unique composition matches closed form")
check(abs(fit2$common) < 1e-12, "orthogonal blocks give exactly zero common variance")
check(abs(fit2$r2_full - (a^2 + b^2) / tss) < 1e-12, "full R-squared matches closed form")

# ------------------------------------------------------------- 3. suppression
# Two strongly correlated predictors with opposite-signed effects. The
# commonality term must come out NEGATIVE and must not be clamped.
rho <- 0.9
a3 <- rnorm(n)
b3 <- rho * a3 + sqrt(1 - rho^2) * rnorm(n)
y3 <- a3 - b3 + 0.1 * rnorm(n)
fit3 <- t5_commonality_cohort(
  matrix(y3, ncol = 1L), matrix(1, nrow = n, ncol = 1L),
  matrix(a3, ncol = 1L), matrix(b3, ncol = 1L)
)
check(fit3$common < -0.5, sprintf("suppression yields a negative common term (%.3f)", fit3$common))
check(abs(fit3$r2_covariates + fit3$unique_fibrosis + fit3$unique_composition +
            fit3$common - fit3$r2_full) < 1e-12,
      "additivity holds under suppression")

# ------------------------------------------------------ 4. agreement with lm()
d <- data.frame(y = Y[, 1L], sex = C[, "sex"], nas = C[, "nas"], fib = A[, 1L],
                pc1 = B[, 1L], pc2 = B[, 2L], pc3 = B[, 3L])
full <- stats::lm(y ~ sex + nas + fib + pc1 + pc2 + pc3, data = d)
reduced <- stats::lm(y ~ sex + nas + pc1 + pc2 + pc3, data = d)
check(abs(summary(full)$r.squared - fit$r2_full[[1L]]) < 1e-12, "full R-squared matches lm()")
anova_fib <- stats::anova(reduced, full)
check(abs(anova_fib$F[2L] - fit$f_unique_fibrosis[[1L]]) < 1e-9,
      "partial F for fibrosis matches anova(lm)")
reduced_b <- stats::lm(y ~ sex + nas + fib, data = d)
check(abs(stats::anova(reduced_b, full)$F[2L] - fit$f_unique_composition[[1L]]) < 1e-9,
      "partial F for the composition block matches anova(lm)")

# ------------------------------------------------- 5. precomputed-parts branch
parts <- t5_baseline_parts(Y, C, A)
fit_fast <- t5_commonality_cohort(Y, C, A, B, parts = parts)
check(identical(fit_fast, fit), "the precomputed-parts fast path is byte-identical to the slow path")

# ---------------------------------------------------------- 5. rank deficiency
Bdup <- cbind(B, B[, 1L, drop = FALSE])
check(is.null(t5_rss(Y, cbind(C, A, Bdup))), "rank-deficient design refuses rather than projecting")

# ------------------------------------------------------------- 6. aggregation
# Aggregation is a linear map, so the identity must survive it.
per_cohort <- list(
  t5_commonality_cohort(Y[1:150, , drop = FALSE], C[1:150, ], A[1:150, , drop = FALSE],
                        B[1:150, , drop = FALSE]),
  t5_commonality_cohort(Y[151:300, , drop = FALSE], C[151:300, ], A[151:300, , drop = FALSE],
                        B[151:300, , drop = FALSE])
)
agg <- t5_aggregate_cohorts(per_cohort)
resid_agg <- with(agg, r2_covariates + unique_fibrosis + unique_composition + common - r2_full)
check(max(abs(resid_agg)) < 1e-12,
      sprintf("aggregated components sum to the aggregated full R-squared (max |error| %.3g)",
              max(abs(resid_agg))))
check(abs(agg$n_donors - 300) < 1e-12, "aggregation counts every donor once")

# ---------------------------------------------------- 7. fibrosis increment id
check(max(abs(agg$fibrosis_incremental - (agg$unique_fibrosis + agg$common))) < 1e-12,
      "fibrosis increment equals unique fibrosis plus common")

# ------------------------------------------------------- 8. coupling fraction
frac <- t5_coupling_fraction(agg$common, agg$fibrosis_incremental)
check(max(abs(frac - agg$common / agg$fibrosis_incremental)) < 1e-12,
      "raw coupling fraction is common over the fibrosis increment")
corrected <- t5_coupling_fraction(agg$common, agg$fibrosis_incremental,
                                  common_null_mean = 0.01)
check(all(corrected < frac), "null-mean correction lowers the coupling fraction")
check(any(t5_coupling_fraction(c(-0.02, 0.02), c(0.1, 0.1)) < 0),
      "a negative common term produces a negative fraction rather than a clamped zero")

# ---------------------------------------------------------- 9. block permuting
group <- rep(c("c1", "c2"), each = 150L)
permuted <- t5_permute_block_within(B, group)
check(!identical(permuted, B), "permutation actually moves rows")
for (g in unique(group)) {
  rows <- which(group == g)
  check(identical(sort(permuted[rows, 1L]), sort(B[rows, 1L])),
        paste0("permutation stays inside cohort ", g))
  check(identical(order(permuted[rows, 1L]), order(permuted[rows, 1L])) &&
          all(vapply(seq_len(nrow(permuted[rows, , drop = FALSE])), function(i) {
            j <- which(abs(B[rows, 1L] - permuted[rows, , drop = FALSE][i, 1L]) < 1e-12)
            length(j) >= 1L && any(abs(B[rows, , drop = FALSE][j, 2L] -
                                         permuted[rows, , drop = FALSE][i, 2L]) < 1e-12)
          }, logical(1))),
        paste0("permutation keeps a donor's composition columns together in ", g))
}
strata <- rep(rep(c("s1", "s2", "s3"), each = 50L), 2L)
permuted_strat <- t5_permute_block_within(B, group, strata)
key <- paste(group, strata)
check(all(vapply(unique(key), function(k) {
  rows <- which(key == k)
  identical(sort(permuted_strat[rows, 1L]), sort(B[rows, 1L]))
}, logical(1))), "stratified permutation stays inside cohort-by-stratum cells")

message("\nT5 TESTS PASSED: ", passed, " checks")
