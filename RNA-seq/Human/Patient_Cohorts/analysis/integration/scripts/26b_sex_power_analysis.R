#!/usr/bin/env Rscript
# 26b_sex_power_analysis.R
# ---------------------------------------------------------------------------
# Formal power analysis for the sex-stratified DE interaction model.
#
# Reviewer concern: the interaction test (101 F controls, 52 M controls)
# may be underpowered, and XIST (used for sex inference) may change with
# liver disease, biasing inference in unannotated cohorts.
#
# Three analyses:
#   Task 1: Formal power analysis via simulation (no pwr package needed)
#   Task 2: XIST disease-dependence check (annotated-sex cohorts only)
#   Task 3: k-means sex inference robustness (silhouette, per-cohort)
#
# Output:
#   RNA-seq/results/audit_sensitivity/sex_power_analysis/
#     REPORT.md, power_curves.csv, xist_disease_correlation.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning ----
set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(yaml)
  library(cluster)  # silhouette
  library(stats)
})

cat("======================================================\n")
cat("26b_sex_power_analysis.R — Formal sex-interaction power analysis\n")
cat("======================================================\n\n")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE     <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT      <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR     <- file.path(INT, "results/integration")
OUT      <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_power_analysis")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Gene IDs (GENCODE v49, versioned Ensembl)
XIST_ID  <- "ENSG00000229807.15"
DDX3Y_ID <- "ENSG00000012817.17"

# ============================================================================
# Load data
# ============================================================================
cat("===== Loading data =====\n")
dge  <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc   <- fread(file.path(INT, "qc/sample_qc_report.csv"))

# Mega cohorts
ycfg <- read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("Mega cohorts:", paste(mega, collapse = ", "), "\n")

# Filter to mega + pass_technical + pass_sex (same as Script 26)
pass_ids <- qc[pass_technical == TRUE & pass_sex == TRUE, sample_id]
keep <- colnames(dge) %in% pass_ids & dge$samples$dataset %in% mega
dge_mega <- dge[, keep]

# Attach inferred sex
matched_sex <- meta$inferred_sex[match(colnames(dge_mega), meta$sample_id)]
dge_mega$samples$inferred_sex <- matched_sex

# Drop NA sex
na_sex <- is.na(dge_mega$samples$inferred_sex)
if (any(na_sex)) {
  cat("Dropping", sum(na_sex), "samples with NA sex\n")
  dge_mega <- dge_mega[, !na_sex]
}

cat("Final sample count:", ncol(dge_mega), "\n\n")

# Group counts
info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  inferred_sex = factor(dge_mega$samples$inferred_sex),
  dataset      = factor(dge_mega$samples$dataset)
)
rownames(info) <- colnames(dge_mega)

cat("Cell counts for interaction model:\n")
cell_counts <- table(info$group_binary, info$inferred_sex)
print(cell_counts)
cat("\n")

n_fc <- cell_counts["Control", "F"]  # 101
n_mc <- cell_counts["Control", "M"]  # 52
n_fd <- cell_counts["Disease", "F"]  # ~371
n_md <- cell_counts["Disease", "M"]  # ~309
n_total <- sum(cell_counts)
cat(sprintf("F-Control=%d, M-Control=%d, F-Disease=%d, M-Disease=%d, Total=%d\n",
            n_fc, n_mc, n_fd, n_md, n_total))

# ============================================================================
# TASK 1: FORMAL POWER ANALYSIS (simulation-based)
# ============================================================================
cat("\n======================================================\n")
cat("TASK 1: Formal power analysis for sex-by-disease interaction\n")
cat("======================================================\n\n")

# ---------------------------------------------------------------------------
# Approach: Analytic power for a 2x2 interaction in a linear model.
#
# The interaction coefficient beta_int tests:
#   (Disease_F - Control_F) - (Disease_M - Control_M)
#
# Under a fixed-effects linear model, the variance of beta_int is:
#   Var(beta_int) = sigma^2 * (1/n_fc + 1/n_mc + 1/n_fd + 1/n_md)
#
# The test statistic t = beta_int / SE(beta_int) follows t(df_residual).
# For the dream mixed model with random intercepts for dataset, the
# effective residual variance is reduced by the ICC, but we compute the
# conservative fixed-effect power as a LOWER BOUND.
#
# BH correction at m=34,453 genes: effective alpha ~ 0.05/m * pi_0_hat.
# We use the Bonferroni bound (most conservative) and also the BH-aware
# power formula.
# ---------------------------------------------------------------------------

n_genes <- nrow(dge)  # dynamic: 39,276 under kallisto, was 34,453 under STAR
alpha_nom <- 0.05

# Effective sample size factor for interaction
# SE(beta_int)^2 = sigma^2 * (1/n1 + 1/n2 + 1/n3 + 1/n4)
inv_ess <- 1/n_fc + 1/n_mc + 1/n_fd + 1/n_md
ess_int <- 1 / inv_ess  # effective sample size for interaction
cat(sprintf("Effective sample size (interaction): %.1f\n", ess_int))
cat(sprintf("  1/n_fc + 1/n_mc + 1/n_fd + 1/n_md = %.5f\n", inv_ess))

# Residual df (fixed-effects model: N - 4 cells - (k-1) datasets)
n_datasets <- nlevels(info$dataset)
df_residual <- n_total - 4 - (n_datasets - 1)
cat(sprintf("Residual df (fixed model): %d\n", df_residual))

# For BH correction, use Bonferroni-equivalent threshold as conservative bound
# and also the BH-empirical threshold
alpha_bonf <- alpha_nom / n_genes
cat(sprintf("Bonferroni threshold: %.2e\n", alpha_bonf))

# BH threshold depends on the rank of the p-value. For power calculation,
# we use the "average" BH threshold: if pi_0 fraction of genes are null
# and we want padj < 0.05, the effective per-test alpha is ~ 0.05 * pi_0.
# Conservative: assume pi_0 = 1 (all null), so BH ~ Bonferroni.
# Realistic: use pi_0 from existing interaction results.
int_res <- fread(file.path(RDIR, "sex_interaction_dream.csv"))
pi_0_hat <- sum(int_res$padj > 0.5, na.rm = TRUE) / nrow(int_res)
cat(sprintf("Estimated pi_0 (fraction truly null): %.3f\n", pi_0_hat))

# BH effective alpha for a non-null gene ranked near the median significant gene
# Under BH, the k-th smallest p-value is compared to k*alpha/m.
# For power at a "typical" non-null gene: alpha_eff ~ alpha * k_sig / m
# Conservative: alpha_eff = alpha / m (Bonferroni)
# We report both Bonferroni and BH-realistic

# ---------------------------------------------------------------------------
# Power curves: sweep Cohen's d from 0.1 to 1.5
# ---------------------------------------------------------------------------
d_grid <- seq(0.05, 1.5, by = 0.05)

# Convert Cohen's d to the interaction logFC scale:
# d = beta_int / sigma_residual
# So beta_int = d * sigma_residual, and
# t = beta_int / SE = d * sigma / sqrt(sigma^2 * inv_ess) = d / sqrt(inv_ess)
# = d * sqrt(ess_int)
#
# This is the noncentrality parameter.

power_results <- data.table(
  cohens_d = d_grid,
  ncp = d_grid * sqrt(ess_int),
  power_bonf = NA_real_,
  power_bh_realistic = NA_real_
)

# BH-realistic: assume ~5% of genes are truly differentially expressed
# (conservative for sex interaction). Then alpha_eff ~ 0.05 * 0.05 / 0.95 ~= 0.05/m * correction
# Actually, under BH with pi_0, the effective per-comparison alpha for the
# weakest-rejected gene is approximately alpha * (1 - pi_0) / pi_0 + alpha/m.
# For simplicity: use alpha_bh = alpha * n_sig_expected / m
# But the cleaner way: the critical t for BH at a GIVEN gene is the t that
# gives p = padj_target * rank / m. For the MEDIAN non-null gene:
# If n_true = 1000, median rank ~ 500, then p_crit = 0.05 * 500 / 34453.
# We'll sweep multiple scenarios.

# Simple approach: compute power at Bonferroni threshold and at a BH-equivalent
# threshold assuming 1000 true positives (optimistic) and 100 (pessimistic)
for (i in seq_along(d_grid)) {
  ncp_i <- power_results$ncp[i]

  # t critical value for Bonferroni
  t_crit_bonf <- qt(1 - alpha_bonf / 2, df = df_residual)
  power_results$power_bonf[i] <- 1 - pt(t_crit_bonf, df = df_residual, ncp = ncp_i) +
                                      pt(-t_crit_bonf, df = df_residual, ncp = ncp_i)

  # BH-realistic: median non-null gene at rank ~ 500 out of 34453
  # p_crit = 0.05 * 500 / 34453 = 7.26e-4
  alpha_bh_500 <- 0.05 * 500 / n_genes
  t_crit_bh <- qt(1 - alpha_bh_500 / 2, df = df_residual)
  power_results$power_bh_realistic[i] <- 1 - pt(t_crit_bh, df = df_residual, ncp = ncp_i) +
                                              pt(-t_crit_bh, df = df_residual, ncp = ncp_i)
}

# Find d where power = 80%
d80_bonf <- approx(power_results$power_bonf, power_results$cohens_d, xout = 0.80)$y
d80_bh   <- approx(power_results$power_bh_realistic, power_results$cohens_d, xout = 0.80)$y

cat("\n===== POWER ANALYSIS RESULTS =====\n")
cat(sprintf("Detectable d at 80%% power (Bonferroni): d >= %.3f\n", d80_bonf))
cat(sprintf("Detectable d at 80%% power (BH k=500):   d >= %.3f\n", d80_bh))

# Power at specific d values
for (d_check in c(0.2, 0.3, 0.5, 0.8)) {
  idx <- which.min(abs(power_results$cohens_d - d_check))
  cat(sprintf("  Power at d=%.1f: Bonferroni=%.1f%%, BH-realistic=%.1f%%\n",
              d_check, 100 * power_results$power_bonf[idx],
              100 * power_results$power_bh_realistic[idx]))
}

# ---------------------------------------------------------------------------
# Simulation-based validation (parametric bootstrap, 5000 reps at key d values)
# ---------------------------------------------------------------------------
cat("\n--- Simulation-based validation ---\n")
n_sim <- 5000
sim_d_vals <- c(0.2, 0.3, 0.5, 0.8, 1.0)
sim_results <- data.table(
  cohens_d = sim_d_vals,
  sim_power_bonf = NA_real_,
  sim_power_bh = NA_real_,
  sim_n_sig_mean = NA_real_,
  sim_n_sig_median = NA_real_
)

cat("Running", n_sim, "simulations per d value...\n")

for (si in seq_along(sim_d_vals)) {
  d_val <- sim_d_vals[si]

  # Simulate p-values: 34,453 genes, 1 is truly non-null with effect d
  # Under H0: p ~ Uniform(0,1)
  # Under H1: p from non-central t with ncp = d * sqrt(ess_int)
  ncp_val <- d_val * sqrt(ess_int)

  bonf_hits <- 0
  bh_hits <- 0
  n_sig_vec <- numeric(n_sim)

  for (r in 1:n_sim) {
    # Generate t-statistics for n_genes null genes + the focal gene
    t_null <- rt(n_genes - 1, df = df_residual)
    t_alt  <- rt(1, df = df_residual, ncp = ncp_val)

    # p-values (two-sided)
    p_null <- 2 * pt(-abs(t_null), df = df_residual)
    p_alt  <- 2 * pt(-abs(t_alt),  df = df_residual)

    all_p <- c(p_null, p_alt)
    all_padj <- p.adjust(all_p, method = "BH")

    # Did the alternative gene pass Bonferroni?
    bonf_hits <- bonf_hits + (p_alt < alpha_bonf)

    # Did the alternative gene pass BH?
    bh_hits <- bh_hits + (all_padj[n_genes] < 0.05)

    # How many total genes passed BH?
    n_sig_vec[r] <- sum(all_padj < 0.05)
  }

  sim_results$sim_power_bonf[si]   <- bonf_hits / n_sim
  sim_results$sim_power_bh[si]     <- bh_hits / n_sim
  sim_results$sim_n_sig_mean[si]   <- mean(n_sig_vec)
  sim_results$sim_n_sig_median[si] <- median(n_sig_vec)

  cat(sprintf("  d=%.1f: sim_power_bonf=%.1f%%, sim_power_bh=%.1f%%, mean_n_sig=%.1f\n",
              d_val, 100 * sim_results$sim_power_bonf[si],
              100 * sim_results$sim_power_bh[si],
              sim_results$sim_n_sig_mean[si]))
}

# ---------------------------------------------------------------------------
# More realistic simulation: ~2000 genes have true effects (matching the
# observed ~2653 padj<0.05 genes from the interaction model)
# ---------------------------------------------------------------------------
cat("\n--- Realistic multi-gene simulation (2000 true effects) ---\n")
n_true <- 2000
n_null_sim <- n_genes - n_true

realistic_sim <- data.table(
  cohens_d_true = c(0.2, 0.3, 0.4, 0.5),
  avg_power_per_gene = NA_real_,
  n_detected_mean = NA_real_,
  n_detected_sd = NA_real_,
  fdr_mean = NA_real_
)

n_sim_realistic <- 1000

for (si in seq_along(realistic_sim$cohens_d_true)) {
  d_val <- realistic_sim$cohens_d_true[si]
  ncp_val <- d_val * sqrt(ess_int)

  detected_vec <- numeric(n_sim_realistic)
  fdr_vec <- numeric(n_sim_realistic)
  tp_vec <- numeric(n_sim_realistic)

  for (r in 1:n_sim_realistic) {
    t_null <- rt(n_null_sim, df = df_residual)
    t_alt  <- rt(n_true,     df = df_residual, ncp = ncp_val)

    p_null <- 2 * pt(-abs(t_null), df = df_residual)
    p_alt  <- 2 * pt(-abs(t_alt),  df = df_residual)

    all_p <- c(p_null, p_alt)
    all_padj <- p.adjust(all_p, method = "BH")

    is_sig <- all_padj < 0.05
    n_sig_total <- sum(is_sig)
    n_tp <- sum(is_sig[(n_null_sim + 1):n_genes])  # true positives
    n_fp <- sum(is_sig[1:n_null_sim])               # false positives

    detected_vec[r] <- n_sig_total
    tp_vec[r] <- n_tp
    fdr_vec[r] <- if (n_sig_total > 0) n_fp / n_sig_total else 0
  }

  realistic_sim$avg_power_per_gene[si] <- mean(tp_vec) / n_true
  realistic_sim$n_detected_mean[si]    <- mean(detected_vec)
  realistic_sim$n_detected_sd[si]      <- sd(detected_vec)
  realistic_sim$fdr_mean[si]           <- mean(fdr_vec)

  cat(sprintf("  d=%.1f: avg_power=%.1f%%, n_detected=%.0f+/-%.0f, FDR=%.3f\n",
              d_val,
              100 * realistic_sim$avg_power_per_gene[si],
              realistic_sim$n_detected_mean[si],
              realistic_sim$n_detected_sd[si],
              realistic_sim$fdr_mean[si]))
}

# Save power curves
fwrite(power_results, file.path(OUT, "power_curves.csv"))
fwrite(sim_results, file.path(OUT, "simulation_power_validation.csv"))
fwrite(realistic_sim, file.path(OUT, "realistic_multi_gene_simulation.csv"))
cat("\nSaved: power_curves.csv, simulation_power_validation.csv, realistic_multi_gene_simulation.csv\n")


# ============================================================================
# TASK 2: XIST DISEASE-DEPENDENCE CHECK
# ============================================================================
cat("\n======================================================\n")
cat("TASK 2: XIST expression versus disease stage (annotated cohorts)\n")
cat("======================================================\n\n")

# Work with the FULL dataset (not just mega) to maximize annotated-sex cohorts
# Annotated-sex cohorts: GSE126848, GSE130970, GSE162694, GSE174478, GSE167523,
#                        GSE193066, PRJNA512027
# Inferred-sex cohorts: GSE135251, GSE213621, GSE240729

# Compute CPM for XIST and DDX3Y
if (XIST_ID %in% rownames(dge) && DDX3Y_ID %in% rownames(dge)) {
  cat("Computing CPM for XIST and DDX3Y across all samples...\n")
  cpm_mat <- cpm(dge, log = FALSE)
  xist_cpm  <- cpm_mat[XIST_ID, ]
  ddx3y_cpm <- cpm_mat[DDX3Y_ID, ]

  sex_df <- data.table(
    sample_id    = colnames(dge),
    dataset      = dge$samples$dataset,
    group_binary = dge$samples$group_binary,
    xist_cpm     = xist_cpm,
    xist_log1p   = log1p(xist_cpm),
    ddx3y_cpm    = ddx3y_cpm,
    ddx3y_log1p  = log1p(ddx3y_cpm)
  )

  # Merge with meta for annotated sex and disease details
  meta_dt <- as.data.table(meta)
  sex_df <- merge(sex_df, meta_dt[, .(sample_id, sex, inferred_sex, sex_source,
                                       fibrosis_stage, diagnosis_harmonized)],
                  by = "sample_id", all.x = TRUE)

  # --- 2a. XIST vs disease group in ANNOTATED-sex cohorts ---
  annotated <- sex_df[sex_source == "annotated"]
  cat("Annotated-sex samples:", nrow(annotated), "\n")
  cat("Datasets:", paste(unique(annotated$dataset), collapse = ", "), "\n\n")

  # Focus on females (XIST expressors)
  annotated_f <- annotated[sex == "F"]
  annotated_m <- annotated[sex == "M"]

  cat("--- XIST in annotated females by disease group ---\n")
  cat(sprintf("  Female Control (n=%d): XIST log1p-CPM = %.2f +/- %.2f\n",
              sum(annotated_f$group_binary == "Control"),
              mean(annotated_f[group_binary == "Control", xist_log1p]),
              sd(annotated_f[group_binary == "Control", xist_log1p])))
  cat(sprintf("  Female Disease (n=%d): XIST log1p-CPM = %.2f +/- %.2f\n",
              sum(annotated_f$group_binary == "Disease"),
              mean(annotated_f[group_binary == "Disease", xist_log1p]),
              sd(annotated_f[group_binary == "Disease", xist_log1p])))

  # Wilcoxon test: XIST in female controls vs female disease
  wt_f <- wilcox.test(annotated_f[group_binary == "Control", xist_log1p],
                      annotated_f[group_binary == "Disease", xist_log1p])
  cat(sprintf("  Wilcoxon p = %.3e\n", wt_f$p.value))

  cat("\n--- XIST in annotated males by disease group ---\n")
  cat(sprintf("  Male Control (n=%d): XIST log1p-CPM = %.2f +/- %.2f\n",
              sum(annotated_m$group_binary == "Control"),
              mean(annotated_m[group_binary == "Control", xist_log1p]),
              sd(annotated_m[group_binary == "Control", xist_log1p])))
  cat(sprintf("  Male Disease (n=%d): XIST log1p-CPM = %.2f +/- %.2f\n",
              sum(annotated_m$group_binary == "Disease"),
              mean(annotated_m[group_binary == "Disease", xist_log1p]),
              sd(annotated_m[group_binary == "Disease", xist_log1p])))

  wt_m <- wilcox.test(annotated_m[group_binary == "Control", xist_log1p],
                      annotated_m[group_binary == "Disease", xist_log1p])
  cat(sprintf("  Wilcoxon p = %.3e\n", wt_m$p.value))

  # --- 2b. Per-cohort XIST vs disease stage (where fibrosis_stage available) ---
  cat("\n--- Per-cohort XIST correlation with fibrosis stage ---\n")
  xist_cors <- list()

  for (ds in sort(unique(annotated$dataset))) {
    sub_f <- annotated[dataset == ds & sex == "F" & !is.na(fibrosis_stage)]
    if (nrow(sub_f) >= 10 && length(unique(sub_f$fibrosis_stage)) >= 2) {
      # Convert fibrosis stage to numeric
      fstage <- suppressWarnings(as.numeric(gsub("F", "", sub_f$fibrosis_stage)))
      if (sum(!is.na(fstage)) >= 10) {
        cor_test <- cor.test(fstage, sub_f$xist_log1p, method = "spearman",
                             exact = FALSE)
        xist_cors[[ds]] <- data.table(
          dataset = ds,
          sex = "F",
          n = sum(!is.na(fstage)),
          rho = cor_test$estimate,
          p_value = cor_test$p.value,
          xist_mean = mean(sub_f$xist_log1p),
          xist_sd = sd(sub_f$xist_log1p)
        )
        cat(sprintf("  %s (F, n=%d): rho=%.3f, p=%.3e\n",
                    ds, sum(!is.na(fstage)), cor_test$estimate, cor_test$p.value))
      }
    }
    # Also males
    sub_m <- annotated[dataset == ds & sex == "M" & !is.na(fibrosis_stage)]
    if (nrow(sub_m) >= 10 && length(unique(sub_m$fibrosis_stage)) >= 2) {
      fstage_m <- suppressWarnings(as.numeric(gsub("F", "", sub_m$fibrosis_stage)))
      if (sum(!is.na(fstage_m)) >= 10) {
        cor_test_m <- cor.test(fstage_m, sub_m$xist_log1p, method = "spearman",
                               exact = FALSE)
        xist_cors[[paste0(ds, "_M")]] <- data.table(
          dataset = ds,
          sex = "M",
          n = sum(!is.na(fstage_m)),
          rho = cor_test_m$estimate,
          p_value = cor_test_m$p.value,
          xist_mean = mean(sub_m$xist_log1p),
          xist_sd = sd(sub_m$xist_log1p)
        )
        cat(sprintf("  %s (M, n=%d): rho=%.3f, p=%.3e\n",
                    ds, sum(!is.na(fstage_m)), cor_test_m$estimate,
                    cor_test_m$p.value))
      }
    }
  }

  # Also test using group_binary (available for all samples)
  cat("\n--- Per-cohort XIST vs disease group (Control vs Disease) ---\n")
  xist_group_tests <- list()

  for (ds in sort(unique(annotated$dataset))) {
    for (sx in c("F", "M")) {
      sub <- annotated[dataset == ds & sex == sx]
      n_ctrl <- sum(sub$group_binary == "Control")
      n_dis  <- sum(sub$group_binary == "Disease")
      if (n_ctrl >= 3 && n_dis >= 3) {
        wt <- wilcox.test(sub[group_binary == "Control", xist_log1p],
                          sub[group_binary == "Disease", xist_log1p])
        mean_ctrl <- mean(sub[group_binary == "Control", xist_log1p])
        mean_dis  <- mean(sub[group_binary == "Disease", xist_log1p])
        xist_group_tests[[paste0(ds, "_", sx)]] <- data.table(
          dataset = ds,
          sex = sx,
          n_control = n_ctrl,
          n_disease = n_dis,
          xist_mean_control = mean_ctrl,
          xist_mean_disease = mean_dis,
          xist_delta = mean_dis - mean_ctrl,
          wilcox_p = wt$p.value
        )
        cat(sprintf("  %s (%s): ctrl=%.2f (n=%d), dis=%.2f (n=%d), delta=%.2f, p=%.3e\n",
                    ds, sx, mean_ctrl, n_ctrl, mean_dis, n_dis,
                    mean_dis - mean_ctrl, wt$p.value))
      }
    }
  }

  # Combine and save
  xist_cor_dt <- rbindlist(xist_cors, fill = TRUE)
  xist_group_dt <- rbindlist(xist_group_tests, fill = TRUE)

  xist_out <- rbindlist(list(
    if (nrow(xist_cor_dt) > 0)   cbind(xist_cor_dt, test = "spearman_vs_fstage") else NULL,
    if (nrow(xist_group_dt) > 0) cbind(xist_group_dt, test = "wilcox_ctrl_vs_dis") else NULL
  ), fill = TRUE)
  fwrite(xist_out, file.path(OUT, "xist_disease_correlation.csv"))
  cat("\nSaved: xist_disease_correlation.csv\n")

  # --- 2c. Per-cohort sex inference concordance table ---
  cat("\n--- Sex inference concordance (annotated cohorts) ---\n")
  concordance_table <- annotated[, .(
    n_total = .N,
    n_concordant = sum(sex == inferred_sex),
    n_discordant = sum(sex != inferred_sex),
    concordance_pct = round(100 * sum(sex == inferred_sex) / .N, 2)
  ), by = dataset]
  concordance_table <- concordance_table[order(dataset)]
  print(concordance_table)
  fwrite(concordance_table, file.path(OUT, "sex_inference_concordance.csv"))
  cat("Saved: sex_inference_concordance.csv\n")

} else {
  cat("ERROR: XIST or DDX3Y not found in DGE object\n")
}


# ============================================================================
# TASK 3: K-MEANS SEX INFERENCE ROBUSTNESS
# ============================================================================
cat("\n======================================================\n")
cat("TASK 3: K-means sex inference robustness\n")
cat("======================================================\n\n")

# --- 3a. Silhouette width ---
cat("--- 3a. Silhouette width ---\n")
sex_mat <- cbind(log1p(sex_df$xist_cpm), log1p(sex_df$ddx3y_cpm))
colnames(sex_mat) <- c("XIST_log1p", "DDX3Y_log1p")

# Full-dataset k-means (pooled, same as Script 01)
set.seed(42)
km_pooled <- kmeans(sex_mat, centers = 2, nstart = 25)
cluster_xist_means <- tapply(sex_df$xist_cpm, km_pooled$cluster, mean)
female_cluster <- which.max(cluster_xist_means)
sex_df[, kmeans_pooled := fifelse(km_pooled$cluster == female_cluster, "F", "M")]

# Silhouette
d_mat <- dist(sex_mat)
sil <- silhouette(km_pooled$cluster, d_mat)
sil_summary <- summary(sil)
avg_sil <- sil_summary$avg.width
cat(sprintf("Pooled k-means silhouette width: %.4f\n", avg_sil))
cat(sprintf("  Cluster 1 (n=%d): avg sil = %.4f\n",
            sil_summary$clus.sizes[1], sil_summary$clus.avg.widths[1]))
cat(sprintf("  Cluster 2 (n=%d): avg sil = %.4f\n",
            sil_summary$clus.sizes[2], sil_summary$clus.avg.widths[2]))

# Proportion of samples with negative silhouette (ambiguous)
n_neg_sil <- sum(sil[, "sil_width"] < 0)
cat(sprintf("  Samples with negative silhouette: %d (%.1f%%)\n",
            n_neg_sil, 100 * n_neg_sil / nrow(sex_mat)))

# --- 3b. Per-cohort vs pooled k-means comparison ---
cat("\n--- 3b. Per-cohort vs pooled k-means ---\n")
percohort_results <- list()

for (ds in sort(unique(sex_df$dataset))) {
  sub <- sex_df[dataset == ds]
  if (nrow(sub) < 10) next

  sex_mat_ds <- cbind(log1p(sub$xist_cpm), log1p(sub$ddx3y_cpm))

  set.seed(42)
  km_ds <- kmeans(sex_mat_ds, centers = 2, nstart = 25)
  ds_xist_means <- tapply(sub$xist_cpm, km_ds$cluster, mean)
  ds_female_cluster <- which.max(ds_xist_means)
  sub[, kmeans_percohort := fifelse(km_ds$cluster == ds_female_cluster, "F", "M")]

  # Agreement between per-cohort and pooled
  agree <- sum(sub$kmeans_pooled == sub$kmeans_percohort)
  n_ds <- nrow(sub)

  # Silhouette for per-cohort
  d_ds <- dist(sex_mat_ds)
  sil_ds <- silhouette(km_ds$cluster, d_ds)
  avg_sil_ds <- summary(sil_ds)$avg.width

  percohort_results[[ds]] <- data.table(
    dataset = ds,
    n = n_ds,
    pooled_F = sum(sub$kmeans_pooled == "F"),
    pooled_M = sum(sub$kmeans_pooled == "M"),
    percohort_F = sum(sub$kmeans_percohort == "F"),
    percohort_M = sum(sub$kmeans_percohort == "M"),
    agreement_n = agree,
    agreement_pct = round(100 * agree / n_ds, 2),
    silhouette_pooled = round(avg_sil, 4),
    silhouette_percohort = round(avg_sil_ds, 4)
  )

  # If annotated sex available, compare both to truth
  if (ds %in% unique(annotated$dataset)) {
    ann_sub <- annotated[dataset == ds]
    matched <- merge(sub[, .(sample_id, kmeans_pooled, kmeans_percohort)],
                     ann_sub[, .(sample_id, sex)],
                     by = "sample_id")
    if (nrow(matched) > 0) {
      pool_vs_ann <- sum(matched$kmeans_pooled == matched$sex) / nrow(matched)
      perc_vs_ann <- sum(matched$kmeans_percohort == matched$sex) / nrow(matched)
      percohort_results[[ds]]$pooled_vs_annotated_pct <- round(100 * pool_vs_ann, 2)
      percohort_results[[ds]]$percohort_vs_annotated_pct <- round(100 * perc_vs_ann, 2)
    }
  }

  cat(sprintf("  %s (n=%d): pooled-percohort agreement=%.1f%%, sil_pooled=%.3f, sil_percohort=%.3f\n",
              ds, n_ds, 100 * agree / n_ds, avg_sil, avg_sil_ds))
}

percohort_dt <- rbindlist(percohort_results, fill = TRUE)
fwrite(percohort_dt, file.path(OUT, "kmeans_robustness.csv"))
cat("\nSaved: kmeans_robustness.csv\n")

# --- 3c. nstart sensitivity ---
cat("\n--- 3c. nstart sensitivity ---\n")
nstart_vals <- c(1, 5, 10, 25, 50, 100)
nstart_results <- data.table()

for (ns in nstart_vals) {
  set.seed(42)
  km_ns <- kmeans(sex_mat, centers = 2, nstart = ns)
  ds_xist_means <- tapply(sex_df$xist_cpm, km_ns$cluster, mean)
  ds_female_cluster <- which.max(ds_xist_means)
  sex_ns <- ifelse(km_ns$cluster == ds_female_cluster, "F", "M")

  agree_25 <- sum(sex_ns == sex_df$kmeans_pooled)
  nstart_results <- rbind(nstart_results, data.table(
    nstart = ns,
    tot_ss = km_ns$tot.withinss,
    agreement_with_nstart25 = agree_25,
    agreement_pct = round(100 * agree_25 / nrow(sex_df), 4)
  ))
  cat(sprintf("  nstart=%3d: tot.withinss=%.1f, agreement_w_ns25=%.2f%%\n",
              ns, km_ns$tot.withinss, 100 * agree_25 / nrow(sex_df)))
}

fwrite(nstart_results, file.path(OUT, "nstart_sensitivity.csv"))
cat("Saved: nstart_sensitivity.csv\n")


# ============================================================================
# GENERATE REPORT
# ============================================================================
cat("\n======================================================\n")
cat("Generating REPORT.md\n")
cat("======================================================\n\n")

report <- c(
  "# Sex-Stratified DE Power Analysis Report",
  "",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Script: `26b_sex_power_analysis.R`"),
  "",
  "## 1. Study Design",
  "",
  "| Cell | N |",
  "|------|---|",
  sprintf("| Female Control | %d |", n_fc),
  sprintf("| Male Control | %d |", n_mc),
  sprintf("| Female Disease | %d |", n_fd),
  sprintf("| Male Disease | %d |", n_md),
  sprintf("| **Total** | **%d** |", n_total),
  "",
  sprintf("- Genes tested: %d", n_genes),
  sprintf("- Interaction model: `~ group_binary * inferred_sex + (1|dataset)`"),
  sprintf("- Multiple testing: BH at alpha = 0.05"),
  sprintf("- Mega cohorts: %s", paste(mega, collapse = ", ")),
  "",
  "## 2. Power Analysis",
  "",
  "### 2a. Analytic Power (Fixed-Effect Lower Bound)",
  "",
  sprintf("Effective sample size for interaction term: %.1f", ess_int),
  sprintf("Residual df: %d", df_residual),
  "",
  "| Cohen's d | Power (Bonferroni) | Power (BH, k=500) |",
  "|-----------|-------------------|-------------------|"
)

for (i in seq_len(nrow(power_results))) {
  if (power_results$cohens_d[i] %in% c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5)) {
    report <- c(report, sprintf("| %.1f | %.1f%% | %.1f%% |",
                                power_results$cohens_d[i],
                                100 * power_results$power_bonf[i],
                                100 * power_results$power_bh_realistic[i]))
  }
}

report <- c(report,
  "",
  sprintf("**Detectable interaction effect at 80%% power (Bonferroni): d >= %.3f**", d80_bonf),
  sprintf("**Detectable interaction effect at 80%% power (BH k=500): d >= %.3f**", d80_bh),
  "",
  "### 2b. Simulation Validation (5,000 reps, single non-null gene)",
  "",
  "| Cohen's d | Sim Power (Bonferroni) | Sim Power (BH) |",
  "|-----------|----------------------|----------------|"
)

for (i in seq_len(nrow(sim_results))) {
  report <- c(report, sprintf("| %.1f | %.1f%% | %.1f%% |",
                              sim_results$cohens_d[i],
                              100 * sim_results$sim_power_bonf[i],
                              100 * sim_results$sim_power_bh[i]))
}

report <- c(report,
  "",
  "### 2c. Realistic Multi-Gene Simulation (2,000 true effects, 1,000 reps)",
  "",
  "| d (all true) | Avg Power/Gene | N Detected (mean +/- SD) | Mean FDR |",
  "|-------------|---------------|-------------------------|---------|"
)

for (i in seq_len(nrow(realistic_sim))) {
  report <- c(report, sprintf("| %.1f | %.1f%% | %.0f +/- %.0f | %.3f |",
                              realistic_sim$cohens_d_true[i],
                              100 * realistic_sim$avg_power_per_gene[i],
                              realistic_sim$n_detected_mean[i],
                              realistic_sim$n_detected_sd[i],
                              realistic_sim$fdr_mean[i]))
}

report <- c(report,
  "",
  "### 2d. Interpretation",
  "",
  sprintf("The study has %d female and %d male controls. The interaction test", n_fc, n_mc),
  "is powered by the *harmonic mean* of all four cell sizes, making the",
  sprintf("male-control cell (n=%d) the bottleneck.", n_mc),
  "",
  "At d=0.5 (medium effect), per-gene power after multiple testing correction",
  sprintf("is %.1f%% (Bonferroni) to %.1f%% (BH, realistic).",
          100 * power_results$power_bonf[which.min(abs(power_results$cohens_d - 0.5))],
          100 * power_results$power_bh_realistic[which.min(abs(power_results$cohens_d - 0.5))]),
  sprintf("The study can detect interactions of d >= %.2f at 80%% power.", d80_bonf),
  "",
  "## 3. XIST Disease-Dependence",
  ""
)

if (exists("annotated_f") && nrow(annotated_f) > 0) {
  report <- c(report,
    "### 3a. XIST Expression by Disease Group (Annotated Females)",
    "",
    sprintf("- Female Control (n=%d): XIST log1p-CPM = %.2f +/- %.2f",
            sum(annotated_f$group_binary == "Control"),
            mean(annotated_f[group_binary == "Control", xist_log1p]),
            sd(annotated_f[group_binary == "Control", xist_log1p])),
    sprintf("- Female Disease (n=%d): XIST log1p-CPM = %.2f +/- %.2f",
            sum(annotated_f$group_binary == "Disease"),
            mean(annotated_f[group_binary == "Disease", xist_log1p]),
            sd(annotated_f[group_binary == "Disease", xist_log1p])),
    sprintf("- Wilcoxon p = %.3e", wt_f$p.value),
    ""
  )

  if (wt_f$p.value < 0.05) {
    report <- c(report,
      "**WARNING**: XIST expression differs significantly between female controls",
      "and female disease samples. This could bias k-means sex inference in",
      "unannotated cohorts (GSE135251, GSE213621, GSE240729).",
      ""
    )
  } else {
    report <- c(report,
      "XIST expression does NOT differ significantly between female controls",
      "and female disease samples. Sex inference by k-means is unlikely to be",
      "biased by disease status.",
      ""
    )
  }

  report <- c(report,
    "### 3b. Per-Cohort XIST vs Disease Group",
    "",
    "| Dataset | Sex | N Control | N Disease | XIST Delta | Wilcoxon p |",
    "|---------|-----|-----------|-----------|------------|------------|"
  )

  for (i in seq_len(nrow(xist_group_dt))) {
    row <- xist_group_dt[i]
    report <- c(report, sprintf("| %s | %s | %d | %d | %.3f | %.3e |",
                                row$dataset, row$sex, row$n_control, row$n_disease,
                                row$xist_delta, row$wilcox_p))
  }
}

report <- c(report,
  "",
  "### 3c. Per-Cohort Concordance (Annotated vs k-means Inferred)",
  "",
  "| Dataset | N | Concordant | Discordant | Concordance % |",
  "|---------|---|-----------|------------|--------------|"
)

if (exists("concordance_table")) {
  for (i in seq_len(nrow(concordance_table))) {
    row <- concordance_table[i]
    report <- c(report, sprintf("| %s | %d | %d | %d | %.1f%% |",
                                row$dataset, row$n_total, row$n_concordant,
                                row$n_discordant, row$concordance_pct))
  }
}

report <- c(report,
  "",
  "## 4. K-Means Robustness",
  "",
  sprintf("### 4a. Silhouette Width: %.4f", avg_sil),
  "",
  "Interpretation:",
  "- > 0.70: Strong clustering",
  "- 0.50-0.70: Reasonable clustering",
  "- 0.25-0.50: Weak clustering",
  "- < 0.25: No meaningful clustering",
  "",
  sprintf("Samples with negative silhouette (ambiguous sex): %d / %d (%.1f%%)",
          n_neg_sil, nrow(sex_mat), 100 * n_neg_sil / nrow(sex_mat)),
  "",
  "### 4b. Per-Cohort vs Pooled K-Means Agreement",
  "",
  "| Dataset | N | Pooled-PerCohort Agree | Sil (pooled) | Sil (per-cohort) | Pooled vs Ann | PerCoh vs Ann |",
  "|---------|---|----------------------|--------------|-----------------|--------------|--------------|"
)

for (i in seq_len(nrow(percohort_dt))) {
  row <- percohort_dt[i]
  ann_p <- if (!is.null(row$pooled_vs_annotated_pct) && !is.na(row$pooled_vs_annotated_pct)) {
    sprintf("%.1f%%", row$pooled_vs_annotated_pct)
  } else "N/A"
  ann_pc <- if (!is.null(row$percohort_vs_annotated_pct) && !is.na(row$percohort_vs_annotated_pct)) {
    sprintf("%.1f%%", row$percohort_vs_annotated_pct)
  } else "N/A"
  report <- c(report, sprintf("| %s | %d | %.1f%% | %.3f | %.3f | %s | %s |",
                              row$dataset, row$n, row$agreement_pct,
                              row$silhouette_pooled, row$silhouette_percohort,
                              ann_p, ann_pc))
}

report <- c(report,
  "",
  "### 4c. nstart Sensitivity",
  "",
  "| nstart | Total Within-SS | Agreement with nstart=25 |",
  "|--------|----------------|------------------------|"
)

for (i in seq_len(nrow(nstart_results))) {
  row <- nstart_results[i]
  report <- c(report, sprintf("| %d | %.1f | %.2f%% |",
                              row$nstart, row$tot_ss, row$agreement_pct))
}

report <- c(report,
  "",
  "## 5. Conclusion",
  "",
  "This power analysis should be reported in the manuscript supplement as:",
  "",
  sprintf("> The sex-by-disease interaction test was conducted on %d samples", n_total),
  sprintf("> (%d female controls, %d male controls, %d female disease, %d male disease)",
          n_fc, n_mc, n_fd, n_md),
  sprintf("> across %d genes with BH correction. Analytic power calculation shows", n_genes),
  sprintf("> the study is powered to detect interaction effects of Cohen's d >= %.2f", d80_bonf),
  "> at 80% power (Bonferroni threshold). XIST expression was assessed for",
  "> disease-dependence in cohorts with annotated sex to verify sex inference",
  "> reliability.",
  ""
)

writeLines(report, file.path(OUT, "REPORT.md"))
cat("\nSaved: REPORT.md\n")

cat("\n===== DONE =====\n")
cat("All outputs in:", OUT, "\n")
