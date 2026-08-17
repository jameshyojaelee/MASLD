#!/usr/bin/env Rscript
# 09b: the central Q2 adjudication. Is the feature-by-window effect matrix rank 1?
#
# The two hypotheses are structures on M (features x windows), not on individual
# effects:
#
#   H_amp  M = a b^T. One loading vector across features, one scalar gain per
#          window. Every feature moves in the same direction, in fixed proportion,
#          and only the overall magnitude changes along the axis.
#   H_seq  rank >= 2. Features have different shapes, so which programmes are
#          moving depends on where you are along the axis.
#
# The statistic is PVE1, the share of squared singular value carried by the first
# component of the UNCENTERED SVD. Uncentered because H_amp is a b^T with no
# intercept; centering would remove the very structure being tested.
#
# The null is the hard part. A permutation null is wrong here: permuting would
# destroy the rank-1 structure and ask "is M better than noise", which is trivially
# yes. The question is the opposite -- is M WORSE than rank 1 by more than
# measurement error explains. So the null is a parametric bootstrap under H_amp
# itself, using the observed per-entry standard errors. This is the matrix
# generalisation of the disattenuation that produced the r = 0.880/0.965/0.959
# cross-sectional result.

suppressPackageStartupMessages({ library(data.table) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
MATRIX_FILE <- Sys.getenv("CAB_EFFECT_MATRIX", "")
LABEL <- Sys.getenv("CAB_MATRIX_LABEL", "M_gene")
assert_true(nzchar(OUT) && nzchar(MATRIX_FILE), "CAB_OUT_ROOT and CAB_EFFECT_MATRIX are required")
dir.create(file.path(OUT, "q2"), recursive = TRUE, showWarnings = FALSE)

pre <- read_prespec()
N_BOOT <- pre$q2_decision_rule$n_bootstrap
set.seed(pre$seeds$sim_base)

# Long format: feature, window, estimate, se
dt <- fread(MATRIX_FILE)
assert_true(all(c("feature", "window", "estimate", "se") %in% names(dt)),
            "Effect matrix needs columns feature, window, estimate, se")
M <- as.matrix(dcast(dt, feature ~ window, value.var = "estimate")[, -1])
S <- as.matrix(dcast(dt, feature ~ window, value.var = "se")[, -1])
rn <- dcast(dt, feature ~ window, value.var = "estimate")$feature
rownames(M) <- rownames(S) <- rn
ok <- complete.cases(M) & complete.cases(S) & rowSums(S > 0) == ncol(S)
M <- M[ok, , drop = FALSE]; S <- S[ok, , drop = FALSE]
assert_true(nrow(M) >= 10L && ncol(M) >= 3L,
            sprintf("Matrix too small after filtering: %d x %d", nrow(M), ncol(M)))
log_step(LABEL, ": ", nrow(M), " features x ", ncol(M), " windows")

pve1 <- function(X) { d <- svd(X)$d; d[1]^2 / sum(d^2) }

# Weighted rank-1 by alternating least squares, weights 1/se^2. Unweighted SVD
# would let the noisiest entries define the shared structure.
wrank1 <- function(X, W, iters = 200L, tol = 1e-10) {
  s <- svd(X, nu = 1, nv = 1)
  a <- s$u[, 1] * sqrt(s$d[1]); b <- s$v[, 1] * sqrt(s$d[1])
  prev <- Inf
  for (i in seq_len(iters)) {
    for (r in seq_len(nrow(X))) {
      den <- sum(W[r, ] * b^2); a[r] <- if (den > 0) sum(W[r, ] * X[r, ] * b) / den else 0
    }
    for (cc in seq_len(ncol(X))) {
      den <- sum(W[, cc] * a^2); b[cc] <- if (den > 0) sum(W[, cc] * X[, cc] * a) / den else 0
    }
    obj <- sum(W * (X - outer(a, b))^2)
    if (abs(prev - obj) < tol * max(1, abs(prev))) break
    prev <- obj
  }
  list(a = a, b = b, fitted = outer(a, b), obj = obj)
}

W <- 1 / (S^2)
fit <- wrank1(M, W)
# The rank-1 factorisation a b^T is invariant to flipping the signs of BOTH
# vectors, so the raw fit lands on an arbitrary sign. Anchor it once here (gains
# positive) and anchor every bootstrap draw to the same convention below. Without
# this the bootstrap draws straddle both signs and the gain CIs collapse to
# +/-|gain|, which is what the first run reported.
if (sum(fit$b) < 0) { fit$a <- -fit$a; fit$b <- -fit$b }
obs_pve1 <- pve1(M)
obs_pve2 <- { d <- svd(M)$d; d[2]^2 / sum(d^2) }

log_step("parametric bootstrap under H_amp, ", N_BOOT, " draws")
null_pve1 <- numeric(N_BOOT); null_pve2 <- numeric(N_BOOT)
for (b in seq_len(N_BOOT)) {
  Xs <- fit$fitted + matrix(rnorm(length(M), 0, S), nrow(M), ncol(M))
  d <- svd(Xs)$d
  null_pve1[b] <- d[1]^2 / sum(d^2); null_pve2[b] <- d[2]^2 / sum(d^2)
  if (b %% 500 == 0) log_step("  draw ", b, "/", N_BOOT)
}
p5 <- quantile(null_pve1, 0.05)
amp_ok <- obs_pve1 >= p5

# Test 4, gain monotonicity. Under H_amp the fitted b_w IS the gain, so it is read
# straight off the weighted fit rather than re-estimated. A gain that does not
# increase is evidence against clean amplification even if PVE1 is high.
gains <- data.table(window = colnames(M), gain = fit$b)
gains[, gain_rel_last := gain / gain[.N]]
mono <- all(diff(fit$b) > 0) || all(diff(fit$b) < 0)
N_GAIN_BOOT <- min(N_BOOT, 500L)
boot_b <- matrix(NA_real_, N_GAIN_BOOT, ncol(M))
for (b in seq_len(N_GAIN_BOOT)) {
  Xs <- fit$fitted + matrix(rnorm(length(M), 0, S), nrow(M), ncol(M))
  bb <- wrank1(Xs, W, iters = 50L)$b
  # Anchor each draw to the observed convention before taking quantiles.
  if (sum(bb * fit$b) < 0) bb <- -bb
  boot_b[b, ] <- bb
}
gains[, boot_lo := apply(boot_b, 2, quantile, 0.025, na.rm = TRUE)]
gains[, boot_hi := apply(boot_b, 2, quantile, 0.975, na.rm = TRUE)]
# Monotonicity with uncertainty: does each successive gain exceed the previous one
# in the bootstrap, and does the last exceed the first?
gains[, frac_boot_gt_prev := c(NA_real_, vapply(2:ncol(M), function(j)
  mean(boot_b[, j] > boot_b[, j - 1L], na.rm = TRUE), numeric(1)))]
frac_boot_monotone <- mean(apply(boot_b, 1, function(r) all(diff(r) > 0)), na.rm = TRUE)
gain_ratio_last_first <- fit$b[ncol(M)] / fit$b[1]
gain_ratio_boot <- boot_b[, ncol(M)] / boot_b[, 1]

# Per-feature residual: which features, if any, break rank 1. The prespecification
# predicts the EARLIEST window is the sole violator, mirroring F0->F1.
resid <- M - fit$fitted
feat <- data.table(feature = rownames(M),
                   resid_ss = rowSums(resid^2),
                   resid_chisq = rowSums((resid / S)^2),
                   df = ncol(M) - 1L)
feat[, p := pchisq(resid_chisq, df, lower.tail = FALSE)]
feat[, bh := p.adjust(p, "BH")]
setorder(feat, -resid_chisq)
win <- data.table(window = colnames(M),
                  resid_ss = colSums(resid^2),
                  resid_chisq = colSums((resid / S)^2),
                  frac_total_resid_ss = colSums(resid^2) / sum(resid^2))

verdict <- data.table(
  matrix_label = LABEL, n_features = nrow(M), n_windows = ncol(M),
  observed_pve1 = obs_pve1, observed_pve2 = obs_pve2,
  null_pve1_median = median(null_pve1), null_pve1_p05 = p5,
  excess_pve1_over_null = obs_pve1 - median(null_pve1),
  excess_pve2_over_null = obs_pve2 - median(null_pve2),
  empirical_p_pve1 = (1 + sum(null_pve1 <= obs_pve1)) / (1 + N_BOOT),
  rank1_sufficient = amp_ok,
  gains_monotone = mono,
  frac_bootstrap_draws_monotone = frac_boot_monotone,
  gain_ratio_last_over_first = gain_ratio_last_first,
  gain_ratio_boot_lo = quantile(gain_ratio_boot, 0.025, na.rm = TRUE),
  gain_ratio_boot_hi = quantile(gain_ratio_boot, 0.975, na.rm = TRUE),
  n_features_breaking_rank1_BH05 = sum(feat$bh < 0.05, na.rm = TRUE),
  frac_features_breaking_rank1 = mean(feat$bh < 0.05, na.rm = TRUE),
  worst_window_by_resid = win$window[which.max(win$resid_chisq)],
  earliest_window_is_worst = identical(win$window[which.max(win$resid_chisq)], win$window[1]),
  n_bootstrap = N_BOOT)

tag <- gsub("[^A-Za-z0-9_]", "_", LABEL)
write_tsv_once(verdict, file.path(OUT, "q2", sprintf("%s_rank1_verdict.tsv", tag)))
write_tsv_once(gains, file.path(OUT, "q2", sprintf("%s_window_gains.tsv", tag)))
write_tsv_once(feat, file.path(OUT, "q2", sprintf("%s_feature_residuals.tsv.gz", tag)))
write_tsv_once(win, file.path(OUT, "q2", sprintf("%s_window_residuals.tsv", tag)))
write_tsv_once(data.table(draw = seq_len(N_BOOT), pve1 = null_pve1, pve2 = null_pve2),
               file.path(OUT, "q2", sprintf("%s_rank1_null.tsv.gz", tag)))
print(verdict)
print(gains)
log_step("Q2_RANK1_COMPLETE ", LABEL, " rank1_sufficient=", amp_ok)
