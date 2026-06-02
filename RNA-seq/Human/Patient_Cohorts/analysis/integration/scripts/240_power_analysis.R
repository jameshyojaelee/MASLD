# 240_power_analysis.R
# Phase 1.6 — Pre-registered power analysis for granular staging pillar.
# Output: RNA-seq/results/granular_staging/power_analysis.csv
# Runs BEFORE F3 sub-state discovery scripts (241–244) commit compute.

suppressPackageStartupMessages({
  library(pROC)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== Phase 1.6 Power Analysis ==\n")
cat("Project root:", PROJECT_ROOT, "\n")
cat("Output dir:", OUT_DIR, "\n\n")

results <- list()

# ---------------------------------------------------------------------------
# Test 1 — F3 sub-state per-cohort discovery (two-proportion split)
# ---------------------------------------------------------------------------
# Discovery cohorts and their F3 sample sizes (from initial exploration).
cohort_n_f3 <- c(GSE213621 = 95, GSE135251 = 54, GSE193066 = 41,
                 GSE174478 = 24, GSE130970 = 14, GSE240729 = 10, GSE162694 = 8)
# (PRJNA512027 (F3=33) was permanently removed from the pipeline 2026-05-15.)

# For each cohort, what is the minimum cluster proportion split detectable
# at α=0.05, 1-β=0.80, given two-proportion test under null p1=p2=0.5?
# Cohen's h = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)).
# We solve: given n, find smallest |p1-p2| with power >= 0.80.

cohen_h <- function(p1, p2) abs(2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2)))

power_two_prop <- function(n, p1, p2, alpha = 0.05) {
  # Approximation via Cohen's h and z-test on arcsine-transformed proportions.
  h <- cohen_h(p1, p2)
  z_alpha <- qnorm(1 - alpha / 2)
  pwr <- pnorm(h * sqrt(n / 2) - z_alpha) + pnorm(-h * sqrt(n / 2) - z_alpha)
  pwr
}

t1_rows <- list()
for (cohort in names(cohort_n_f3)) {
  n <- cohort_n_f3[cohort]
  # Search smallest detectable split p2 (with p1 = 1 - p2, p2 < 0.5).
  splits <- seq(0.50, 0.20, by = -0.01)
  pwrs <- sapply(splits, function(p2) power_two_prop(n, 1 - p2, p2))
  idx <- which(pwrs >= 0.80)[1]
  min_split <- if (!is.na(idx)) splits[idx] else NA_real_
  t1_rows[[cohort]] <- data.frame(
    test = "f3_substate_per_cohort_split",
    cohort = cohort,
    n_f3 = n,
    detection_threshold = paste0(round((1 - min_split) * 100), "/", round(min_split * 100)),
    detection_value = if (is.na(min_split)) NA else 1 - min_split,
    power_achieved = if (is.na(idx)) max(pwrs) else pwrs[idx],
    notes = "Min imbalanced split (p1/p2) detectable at α=0.05, 1-β=0.80",
    stringsAsFactors = FALSE
  )
}
results$test1_per_cohort_split <- do.call(rbind, t1_rows)

cat("Test 1 — Per-cohort cluster split detection:\n")
print(results$test1_per_cohort_split[, c("cohort", "n_f3", "detection_threshold", "power_achieved")], row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# Test 2 — Bimodality detection (Hartigan-style dip test analog via simulation)
# ---------------------------------------------------------------------------
# Simulate two-component Gaussian mixture with separation Δμ/σ ranging.
# Compute Bayesian Information Criterion (BIC) ratio for k=2 vs k=1.
# A simple manual implementation since mclust is unavailable.

fit_gmm_2 <- function(x, max_iter = 100, tol = 1e-6) {
  # 2-component Gaussian mixture via EM. Returns list with mu1, mu2, sd1, sd2, p1, ll.
  n <- length(x)
  # Init by sorted halves
  ord <- sort(x)
  mu1 <- mean(ord[1:floor(n/2)])
  mu2 <- mean(ord[(floor(n/2)+1):n])
  sd1 <- sd(ord[1:floor(n/2)])
  sd2 <- sd(ord[(floor(n/2)+1):n])
  p1 <- 0.5
  ll_old <- -Inf
  for (iter in 1:max_iter) {
    # E-step
    d1 <- p1 * dnorm(x, mu1, sd1)
    d2 <- (1 - p1) * dnorm(x, mu2, sd2)
    g1 <- d1 / (d1 + d2 + 1e-300)
    # M-step
    p1 <- mean(g1)
    mu1 <- sum(g1 * x) / sum(g1)
    mu2 <- sum((1 - g1) * x) / sum(1 - g1)
    sd1 <- sqrt(sum(g1 * (x - mu1)^2) / sum(g1))
    sd2 <- sqrt(sum((1 - g1) * (x - mu2)^2) / sum(1 - g1))
    sd1 <- max(sd1, 1e-6); sd2 <- max(sd2, 1e-6)
    ll <- sum(log(p1 * dnorm(x, mu1, sd1) + (1 - p1) * dnorm(x, mu2, sd2) + 1e-300))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, p1 = p1, ll = ll)
}

bic_diff <- function(x) {
  n <- length(x)
  ll1 <- sum(dnorm(x, mean(x), sd(x), log = TRUE))
  bic1 <- -2 * ll1 + 2 * log(n)  # 2 params: mu, sd
  fit2 <- fit_gmm_2(x)
  bic2 <- -2 * fit2$ll + 5 * log(n)  # 5 params: mu1, mu2, sd1, sd2, p1
  bic1 - bic2  # positive → mixture preferred
}

# Simulate at separations Δμ/σ ∈ {0.5, 1.0, 1.5, 2.0} for n = 279 (pooled F3) and n=95 (largest cohort).
n_sims <- 200
seps <- c(0.5, 1.0, 1.5, 2.0)
ns_test <- c(95, 279)
t2_rows <- list()
for (n in ns_test) {
  for (sep in seps) {
    detect <- 0
    for (i in 1:n_sims) {
      x <- c(rnorm(n / 2, 0, 1), rnorm(n / 2, sep, 1))
      bd <- tryCatch(bic_diff(x), error = function(e) 0)
      if (bd >= 10) detect <- detect + 1
    }
    pwr <- detect / n_sims
    t2_rows[[paste0(n, "_", sep)]] <- data.frame(
      test = "bimodality_detection",
      cohort = ifelse(n == 279, "pooled", "GSE213621_largest"),
      n_f3 = n,
      detection_threshold = sprintf("Δμ/σ=%.1f", sep),
      detection_value = sep,
      power_achieved = pwr,
      notes = "Power to detect 2-Gaussian mixture vs 1-Gaussian via BIC≥10 (n_sim=200)",
      stringsAsFactors = FALSE
    )
  }
}
results$test2_bimodality <- do.call(rbind, t2_rows)
cat("Test 2 — Bimodality detection:\n")
print(results$test2_bimodality[, c("cohort", "n_f3", "detection_threshold", "power_achieved")], row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# Test 3 — LOCO-CV cluster prediction AUROC
# ---------------------------------------------------------------------------
# How well can we distinguish AUROC=X from AUROC=0.5 given LOCO-CV held-out n?
# Use DeLong-style bootstrap on simulated rocs.

power_auroc <- function(n_pos, n_neg, true_auc, n_sims = 500, alpha = 0.05) {
  # Simulate scores under a Gaussian model that gives the desired AUC,
  # bootstrap CI, return power (P CI lower > 0.5).
  detect <- 0
  for (i in 1:n_sims) {
    # Convert AUC to mean separation in normal scores: AUC = pnorm(d/sqrt(2))
    d <- qnorm(true_auc) * sqrt(2)
    pos <- rnorm(n_pos, d, 1)
    neg <- rnorm(n_neg, 0, 1)
    roc_obj <- suppressMessages(pROC::roc(c(rep(1, n_pos), rep(0, n_neg)), c(pos, neg), quiet = TRUE))
    ci <- suppressMessages(pROC::ci.auc(roc_obj, method = "delong"))
    if (ci[1] > 0.5) detect <- detect + 1
  }
  detect / n_sims
}

auc_targets <- c(0.60, 0.65, 0.70, 0.75)
# Held-out cohort sizes (smaller cohorts in LOCO-CV)
loco_ns <- list("GSE174478" = c(pos = 12, neg = 12),
                "GSE130970" = c(pos = 7, neg = 7),
                "GSE240729" = c(pos = 5, neg = 5),
                "GSE162694" = c(pos = 4, neg = 4),
                "pooled_holdout" = c(pos = 30, neg = 30))

t3_rows <- list()
for (held in names(loco_ns)) {
  for (auc in auc_targets) {
    nn <- loco_ns[[held]]
    pwr <- power_auroc(nn["pos"], nn["neg"], auc, n_sims = 200)
    t3_rows[[paste0(held, "_", auc)]] <- data.frame(
      test = "loco_cv_auroc",
      cohort = held,
      n_f3 = sum(nn),
      detection_threshold = sprintf("AUC=%.2f vs 0.5", auc),
      detection_value = auc,
      power_achieved = pwr,
      notes = "Power that DeLong CI lower > 0.5 (n_sim=200)",
      stringsAsFactors = FALSE
    )
  }
}
results$test3_loco_auroc <- do.call(rbind, t3_rows)
cat("Test 3 — LOCO-CV AUROC detection:\n")
print(results$test3_loco_auroc[, c("cohort", "n_f3", "detection_threshold", "power_achieved")], row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# Test 4 — Olink F3 sub-state per-protein t-test
# ---------------------------------------------------------------------------
# 177 staged Olink subjects. Expected F3 = ~30–50. Per F3 sub-state ~15–25.
# 1,461 proteins × Bonferroni → α = 0.05 / 1461 = 3.4e-5.
# What Cohen's d is detectable at 0.80 power?

power_t_test_d <- function(n_per_group, d, alpha) {
  # Two-sample t-test power calculation. Approximation via non-central t.
  df <- 2 * n_per_group - 2
  ncp <- d * sqrt(n_per_group / 2)
  t_crit <- qt(1 - alpha / 2, df)
  power <- 1 - pt(t_crit, df, ncp) + pt(-t_crit, df, ncp)
  power
}

t4_rows <- list()
n_per_group_grid <- c(10, 15, 20, 25, 30)
d_grid <- c(0.5, 0.8, 1.0, 1.5, 2.0)
for (n in n_per_group_grid) {
  for (d in d_grid) {
    pwr_uncorr <- power_t_test_d(n, d, 0.05)
    pwr_bonf <- power_t_test_d(n, d, 0.05 / 1461)
    t4_rows[[paste0(n, "_", d)]] <- data.frame(
      test = "olink_per_protein_t_test",
      cohort = "Olink_F3_subset",
      n_f3 = 2 * n,
      detection_threshold = sprintf("d=%.1f, n_per_grp=%d", d, n),
      detection_value = d,
      power_achieved = pwr_bonf,
      notes = sprintf("Bonferroni α=%.2e (1461 proteins); uncorrected power=%.2f",
                      0.05 / 1461, pwr_uncorr),
      stringsAsFactors = FALSE
    )
  }
}
results$test4_olink <- do.call(rbind, t4_rows)
cat("Test 4 — Olink per-protein t-test (Bonferroni):\n")
print(results$test4_olink[, c("detection_threshold", "power_achieved", "notes")], row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# Test 5 — COLOC PP4>0.5 enrichment per transition (Fisher exact)
# ---------------------------------------------------------------------------
# 364 PP4>0.5 hits in atlas. Per transition, what subset of bulk-DEGs overlap
# with these hits? Need power to detect a Fisher OR ≥ 2 vs background.
# Simulate: total atlas N = 33943, hit_rate = 364/33943 ≈ 0.011.
# Transition gene set sizes range ~500 (small) to ~5000 (large).

simulate_fisher <- function(N, n_hits, n_set, OR, n_sims = 500, alpha = 0.05) {
  # null hit rate
  p_null <- n_hits / N
  # Under OR, P(hit | in set) = OR * p_null / (1 + (OR-1)*p_null)
  p_alt <- OR * p_null / (1 + (OR - 1) * p_null)
  detect <- 0
  for (i in 1:n_sims) {
    set_hits <- rbinom(1, n_set, p_alt)
    nonset_hits <- n_hits - set_hits
    nonset_hits <- max(0, min(N - n_set, nonset_hits))
    set_nonhits <- n_set - set_hits
    nonset_nonhits <- (N - n_set) - nonset_hits
    if (any(c(set_hits, set_nonhits, nonset_hits, nonset_nonhits) < 0)) next
    tab <- matrix(c(set_hits, set_nonhits, nonset_hits, nonset_nonhits), nrow = 2)
    p <- tryCatch(fisher.test(tab, alternative = "greater")$p.value, error = function(e) 1)
    if (p < alpha) detect <- detect + 1
  }
  detect / n_sims
}

t5_rows <- list()
n_set_grid <- c(200, 500, 1000, 2000, 5000)
OR_grid <- c(1.5, 2.0, 3.0, 5.0)
for (n_set in n_set_grid) {
  for (OR in OR_grid) {
    pwr <- simulate_fisher(33943, 364, n_set, OR, n_sims = 200, alpha = 0.05)
    t5_rows[[paste0(n_set, "_", OR)]] <- data.frame(
      test = "coloc_fisher_per_transition",
      cohort = "atlas",
      n_f3 = NA_integer_,
      detection_threshold = sprintf("n_set=%d, OR=%.1f", n_set, OR),
      detection_value = OR,
      power_achieved = pwr,
      notes = "Fisher exact (greater) at α=0.05, 364 hits in 33943",
      stringsAsFactors = FALSE
    )
  }
}
results$test5_coloc_fisher <- do.call(rbind, t5_rows)
cat("Test 5 — COLOC Fisher enrichment:\n")
print(results$test5_coloc_fisher[, c("detection_threshold", "power_achieved")], row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# Test 6 — LINCS reversal per sub-state (effect size on enrichment z)
# ---------------------------------------------------------------------------
# Compound-level reversal scores. Per sub-state, paired comparison of compound
# scores. With 1,107 LINCS compounds, paired Wilcoxon power for shift δ ∈ {0.1, 0.2, 0.5}σ.

power_wilcox_paired <- function(n, delta_sd, n_sims = 500, alpha = 0.05) {
  detect <- 0
  for (i in 1:n_sims) {
    a <- rnorm(n)
    b <- a + rnorm(n, delta_sd, 0.5)  # induced shift with noise
    p <- suppressWarnings(wilcox.test(a, b, paired = TRUE, alternative = "two.sided")$p.value)
    if (p < alpha) detect <- detect + 1
  }
  detect / n_sims
}

t6_rows <- list()
for (delta in c(0.1, 0.2, 0.3, 0.5)) {
  for (alpha_use in c(0.05, 0.05 / 1107)) {
    pwr <- power_wilcox_paired(1107, delta, n_sims = 200, alpha = alpha_use)
    t6_rows[[paste0(delta, "_", alpha_use)]] <- data.frame(
      test = "lincs_reversal_per_substate",
      cohort = "LINCS",
      n_f3 = 1107L,
      detection_threshold = sprintf("Δ=%.1f σ, α=%.2e", delta, alpha_use),
      detection_value = delta,
      power_achieved = pwr,
      notes = "Paired Wilcoxon over 1107 compounds",
      stringsAsFactors = FALSE
    )
  }
}
results$test6_lincs <- do.call(rbind, t6_rows)
cat("Test 6 — LINCS reversal:\n")
print(results$test6_lincs[, c("detection_threshold", "power_achieved")], row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# Aggregate + decision summary
# ---------------------------------------------------------------------------
all_results <- do.call(rbind, results)
out_path <- file.path(OUT_DIR, "power_analysis.csv")
write.csv(all_results, out_path, row.names = FALSE)
cat("Wrote:", out_path, "\n\n")

# Decision summary
cat("== Power-driven decisions (pre-registered) ==\n\n")

cat("D1 — Per-cohort split detection thresholds:\n")
t1_decision <- results$test1_per_cohort_split
for (i in seq_len(nrow(t1_decision))) {
  cat(sprintf("  %-15s n=%2d  → can detect splits more imbalanced than %s (power=%.2f)\n",
              t1_decision$cohort[i], t1_decision$n_f3[i],
              t1_decision$detection_threshold[i], t1_decision$power_achieved[i]))
}
cat("\nD1 RULE: if per-cohort cluster-split power < 0.6 at expected effect size,\n")
cat("    drop that cohort from discovery tier (e.g., GSE162694 n=8 likely too sparse).\n\n")

cat("D2 — Bimodality detection:\n")
t2_decision <- results$test2_bimodality
for (i in seq_len(nrow(t2_decision))) {
  cat(sprintf("  %-20s Δμ/σ=%.1f power=%.2f\n",
              paste0("n=", t2_decision$n_f3[i]),
              t2_decision$detection_value[i],
              t2_decision$power_achieved[i]))
}
cat("\nD2 RULE: pooled n=279 should detect Δμ/σ ≥ 1.0 with power > 0.80.\n")
cat("    If observed effect is Δμ/σ < 1.0 in real data, flag as borderline; consider downgrading.\n\n")

cat("D3 — Olink per-protein t-test:\n")
t4_decision <- subset(results$test4_olink, detection_value %in% c(0.5, 1.0))
for (i in seq_len(nrow(t4_decision))) {
  cat(sprintf("  %s  power=%.3f (Bonf)\n",
              t4_decision$detection_threshold[i], t4_decision$power_achieved[i]))
}
cat("\nD3 RULE: with n_per_group=15 and d=0.8, Bonferroni-corrected power is severely limited.\n")
cat("    Restrict Olink to PRE-SPECIFIED panels (ECM, senescence, inflammation), report effect-size CIs.\n\n")

cat("D4 — COLOC Fisher per transition:\n")
t5_decision <- subset(results$test5_coloc_fisher, detection_value %in% c(2.0, 3.0))
for (i in seq_len(nrow(t5_decision))) {
  cat(sprintf("  %s  power=%.2f\n",
              t5_decision$detection_threshold[i], t5_decision$power_achieved[i]))
}
cat("\nD4 RULE: For OR=2 detection at 0.80 power, need transition gene set ≥ ~1000 genes.\n")
cat("    Smaller transitions need OR≥3 to be detectable.\n\n")

cat("== Power analysis complete. Inspect", out_path, "before running 241–248. ==\n")
