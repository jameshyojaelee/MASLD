#!/usr/bin/env Rscript
# 09c: criterion (ii) of the Q2 decision rule -- change-point dispersion.
#
# Under scalar amplification every feature shares one gain profile, so the window
# where a feature makes its largest jump is the same window for all features
# (whichever gain increment is biggest). Under ordered sequence, features turn on
# at different points, so those change-points spread out.
#
# Measurement error alone disperses change-points, especially for features with
# small effects, so a raw spread is uninterpretable. The null is again the H_amp
# parametric bootstrap: simulate rank-1 structure plus the observed per-entry SEs,
# recompute the spread, and ask whether the observed spread exceeds it.
#
# Because the Q1 gate routed Q2 onto recorded histology bins, "change-point" here is
# the discrete window of largest increment rather than a continuous position. That
# is the correct analogue for five ordinal stages and is stated as such.

suppressPackageStartupMessages({ library(data.table) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
MATRIX_FILE <- Sys.getenv("CAB_EFFECT_MATRIX", "")
LABEL <- Sys.getenv("CAB_MATRIX_LABEL", "M_gene_histology")
assert_true(nzchar(OUT) && nzchar(MATRIX_FILE), "CAB_OUT_ROOT and CAB_EFFECT_MATRIX are required")

pre <- read_prespec()
N_BOOT <- pre$q2_decision_rule$n_bootstrap
set.seed(pre$seeds$sim_base + 7L)

dt <- fread(MATRIX_FILE)
M <- as.matrix(dcast(dt, feature ~ window, value.var = "estimate")[, -1])
S <- as.matrix(dcast(dt, feature ~ window, value.var = "se")[, -1])
rn <- dcast(dt, feature ~ window, value.var = "estimate")$feature
rownames(M) <- rownames(S) <- rn
ok <- complete.cases(M) & complete.cases(S) & rowSums(S > 0) == ncol(S)
M <- M[ok, , drop = FALSE]; S <- S[ok, , drop = FALSE]
log_step(LABEL, ": ", nrow(M), " features x ", ncol(M), " windows")

# Restrict to features with a real effect somewhere. A feature that is flat noise
# has an arbitrary change-point and would inflate the spread under both hypotheses.
z <- abs(M) / S
strong <- rowSums(z > 3) > 0
log_step("features with |z| > 3 in at least one window: ", sum(strong))
assert_true(sum(strong) >= 100L, "Too few features with a detectable effect")

changepoint <- function(X) {
  inc <- cbind(X[, 1], t(apply(X, 1, diff)))   # increment into each window
  max.col(abs(inc), ties.method = "first")
}
spread <- function(cp) c(sd = sd(cp), iqr = IQR(cp), entropy = {
  p <- table(factor(cp, levels = seq_len(ncol(M)))) / length(cp)
  p <- p[p > 0]; -sum(p * log(p))
})

obs_cp <- changepoint(M[strong, , drop = FALSE])
obs <- spread(obs_cp)

wrank1 <- function(X, W, iters = 100L) {
  s <- svd(X, nu = 1, nv = 1)
  a <- s$u[, 1] * sqrt(s$d[1]); b <- s$v[, 1] * sqrt(s$d[1])
  for (i in seq_len(iters)) {
    for (r in seq_len(nrow(X))) { d <- sum(W[r, ] * b^2); a[r] <- if (d > 0) sum(W[r, ] * X[r, ] * b)/d else 0 }
    for (c2 in seq_len(ncol(X))) { d <- sum(W[, c2] * a^2); b[c2] <- if (d > 0) sum(W[, c2] * X[, c2] * a)/d else 0 }
  }
  outer(a, b)
}
Ms <- M[strong, , drop = FALSE]; Ss <- S[strong, , drop = FALSE]
fitted <- wrank1(Ms, 1 / Ss^2)

log_step("H_amp bootstrap, ", N_BOOT, " draws")
null_sd <- numeric(N_BOOT); null_ent <- numeric(N_BOOT)
for (b in seq_len(N_BOOT)) {
  Xs <- fitted + matrix(rnorm(length(Ms), 0, Ss), nrow(Ms), ncol(Ms))
  s <- spread(changepoint(Xs))
  null_sd[b] <- s["sd"]; null_ent[b] <- s["entropy"]
  if (b %% 500 == 0) log_step("  draw ", b, "/", N_BOOT)
}

q95_sd <- quantile(null_sd, 0.95); q95_ent <- quantile(null_ent, 0.95)
tab <- as.data.table(table(factor(obs_cp, levels = seq_len(ncol(M)))))
setnames(tab, c("window_index", "n_features"))
tab[, window := colnames(M)[as.integer(window_index)]]
tab[, frac := n_features / sum(n_features)]

verdict <- data.table(
  matrix_label = LABEL, n_features_tested = sum(strong), n_windows = ncol(M),
  observed_cp_sd = unname(obs["sd"]), null_cp_sd_median = median(null_sd),
  null_cp_sd_q95 = unname(q95_sd),
  observed_cp_entropy = unname(obs["entropy"]),
  null_cp_entropy_median = median(null_ent), null_cp_entropy_q95 = unname(q95_ent),
  empirical_p_sd = (1 + sum(null_sd >= obs["sd"])) / (1 + N_BOOT),
  empirical_p_entropy = (1 + sum(null_ent >= obs["entropy"])) / (1 + N_BOOT),
  criterion_ii_spread_contained = unname(obs["sd"]) <= unname(q95_sd),
  modal_changepoint_window = tab[which.max(n_features), window],
  modal_changepoint_frac = tab[, max(frac)],
  n_bootstrap = N_BOOT)

tag <- gsub("[^A-Za-z0-9_]", "_", LABEL)
write_tsv_once(verdict, file.path(OUT, "q2", sprintf("%s_changepoint_verdict.tsv", tag)))
write_tsv_once(tab, file.path(OUT, "q2", sprintf("%s_changepoint_distribution.tsv", tag)))
write_tsv_once(data.table(draw = seq_len(N_BOOT), sd = null_sd, entropy = null_ent),
               file.path(OUT, "q2", sprintf("%s_changepoint_null.tsv.gz", tag)))
print(verdict); print(tab)
log_step("Q2_CHANGEPOINT_COMPLETE contained=", verdict$criterion_ii_spread_contained)
