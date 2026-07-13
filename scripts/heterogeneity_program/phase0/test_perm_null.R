#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Phase 0b CALIBRATION TEST for perm_null.R — the most important Phase-0 check.
#
# Synthetic cohort-blocked data: NULL genes have equal within-group variance;
# SPIKE genes have inflated variance in group 2 (a true differential-variability
# signal). A correctly-calibrated permutation null must give:
#   (1) NULL genes  → uniform permutation p-values (KS not rejected; FPR≈0.05),
#   (2) SPIKE genes → high power (most reach q<0.10).
# The DV statistic is a vectorised Brown-Forsythe / diffVar-style Levene mean
# difference: z = |x - per(stratum×group) median|, stat = mean(z|grp2)-mean(z|grp1).
# Group labels are permuted WITHIN stratum and z is recomputed each replicate
# (refit-per-replicate — GATE-B), exercising the engine exactly as T1 will.
# Lightweight (toy data) — safe off-node.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/heterogeneity_program/lib/perm_null.R"))

set.seed(1)
S <- 3L; per_grp <- 20L                      # 3 strata × (20 grp1 + 20 grp2) = 120
strata <- factor(rep(paste0("ds", seq_len(S)), each = 2L * per_grp))
group  <- factor(rep(rep(c("g1", "g2"), each = per_grp), times = S))
N <- length(group)
G_null <- 350L; G_spike <- 50L; G <- G_null + G_spike
is_spike <- c(rep(FALSE, G_null), rep(TRUE, G_spike))

# base mean per gene; null = equal sd; spike = grp2 sd inflated ×2.2 (mean unchanged)
mu  <- rnorm(G, 6, 1.5)
sd1 <- runif(G, 0.5, 1.2)                     # group1 sd
sd2 <- ifelse(is_spike, sd1 * 2.2, sd1)       # group2 sd (spike → inflated, SAME mean)
X <- matrix(0, G, N)
g2 <- group == "g2"
for (j in seq_len(N)) X[, j] <- rnorm(G, mu, if (g2[j]) sd2 else sd1)

# vectorised Brown-Forsythe / Levene DV statistic (group2 - group1 mean |dev|)
rowMed <- function(m) if (ncol(m) == 1L) as.numeric(m) else apply(m, 1L, median)
dv_stat <- function(grp) {
  z <- matrix(0, G, N)
  for (s in levels(strata)) for (gg in levels(grp)) {
    idx <- which(strata == s & grp == gg)
    if (length(idx) == 0L) next
    z[, idx] <- abs(X[, idx, drop = FALSE] - rowMed(X[, idx, drop = FALSE]))
  }
  rowMeans(z[, grp == "g2", drop = FALSE]) - rowMeans(z[, grp == "g1", drop = FALSE])
}

cat("Running perm_null (G=", G, " N=", N, " n_perm=500, within-stratum, refit-per-rep)...\n", sep = "")
t0 <- Sys.time()
res <- perm_null(dv_stat, group = group, strata = strata,
                 n_perm = 500L, seed = 42L, alternative = "greater",
                 n_cores = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))))
cat(sprintf("  done in %.1fs\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

p_null  <- res$p_perm[!is_spike]
p_spike <- res$p_perm[is_spike]
q_spike <- res$q_perm[is_spike]

ks <- suppressWarnings(ks.test(p_null, "punif"))
fpr   <- mean(p_null < 0.05)
meanp <- mean(p_null)
power <- mean(q_spike < 0.10)

cat("── NULL genes (n=", G_null, ") — calibration ──\n", sep = "")
cat(sprintf("  mean p_perm = %.3f   (target ≈ 0.50)\n", meanp))
cat(sprintf("  FPR@0.05    = %.3f   (target ≈ 0.05)\n", fpr))
cat(sprintf("  KS vs Uniform p = %.3f  (target > 0.01 = not rejected)\n", ks$p.value))
cat("── SPIKE genes (n=", G_spike, ") — power ──\n", sep = "")
cat(sprintf("  power (q<0.10) = %.3f   (target high)\n", power))

calib_ok <- ks$p.value > 0.01 && fpr <= 0.10 && meanp > 0.40 && meanp < 0.60
power_ok <- power > 0.50
cat("\n=== ", if (calib_ok && power_ok) "PASS" else "FAIL",
    " (calibration ", if (calib_ok) "OK" else "BAD",
    ", power ", if (power_ok) "OK" else "BAD", ") ===\n", sep = "")
if (!(calib_ok && power_ok)) quit(status = 1)
