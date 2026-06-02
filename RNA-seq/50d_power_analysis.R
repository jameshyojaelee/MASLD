#!/usr/bin/env Rscript
# 50d_power_analysis.R
# ---------------------------------------------------------------------------
# COLOC/GWAS Power Analysis (M2)
#
# Pure calculation — no data loading. Computes minimum detectable effect sizes
# at various MAF bins for each GWAS sample size, and estimates how this
# translates to COLOC power (more significant GWAS associations = more power
# for colocalization).
#
# Key question: What effect sizes are detectable at p<5e-8 for:
#   - AFR (N=6,636)
#   - CSA (N=8,876)
#   - BBJ (N=261K)
#   - UKBB (N=361K)
#
# Output: RNA-seq/results/causal_inference/power_analysis/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR    <- Sys.getenv("MASLD_PROJECT_ROOT",
                          "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/power_analysis")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Script 50d: GWAS/COLOC Power Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. GWAS configurations
# ==============================================================================
gwas_configs <- data.table(
  source = c("PanUKBB AFR", "PanUKBB CSA", "BBJ", "UKBB EUR"),
  N = c(6636L, 8876L, 261406L, 361194L),
  ancestry = c("African", "Central/South Asian", "East Asian", "European")
)

# MAF bins
maf_vals <- seq(0.01, 0.50, by = 0.005)

# Significance threshold
alpha <- 5e-8
z_thresh <- qnorm(1 - alpha / 2)  # two-sided

# ==============================================================================
# 2. Minimum detectable beta at each MAF and N
# ==============================================================================
# For a quantitative trait GWAS with N samples:
#   SE(beta) ≈ 1 / sqrt(2 * N * MAF * (1 - MAF))
#   (assuming standardized trait variance = 1)
#   Minimum detectable |beta| at power ~50% ≈ z_thresh * SE
#   For 80% power, multiply by ~1.4 (adds z_0.8 = 0.84 to z_thresh)

cat("--- Computing minimum detectable betas ---\n")

results <- rbindlist(lapply(seq_len(nrow(gwas_configs)), function(i) {
  N <- gwas_configs$N[i]
  src <- gwas_configs$source[i]

  data.table(
    source = src,
    N = N,
    MAF = maf_vals,
    SE = 1 / sqrt(2 * N * maf_vals * (1 - maf_vals)),
    min_beta_50pct = z_thresh / sqrt(2 * N * maf_vals * (1 - maf_vals)),
    min_beta_80pct = (z_thresh + qnorm(0.80)) / sqrt(2 * N * maf_vals * (1 - maf_vals))
  )
}))

# Variance explained (R^2) at min detectable beta
results[, R2_50pct := 2 * MAF * (1 - MAF) * min_beta_50pct^2]
results[, R2_80pct := 2 * MAF * (1 - MAF) * min_beta_80pct^2]

cat("  Computed power curves for", nrow(gwas_configs), "GWAS sources\n")

# ==============================================================================
# 3. Summary table at key MAF bins
# ==============================================================================
cat("\n--- Summary at key MAF bins ---\n")
key_mafs <- c(0.02, 0.05, 0.10, 0.20, 0.50)
summary_tab <- results[MAF %in% key_mafs,
                        .(source, N, MAF,
                          min_beta_80pct = round(min_beta_80pct, 4),
                          R2_80pct = round(R2_80pct, 6))]
setorder(summary_tab, MAF, -N)
print(summary_tab)

# Interpretation
cat("\n  Interpretation:\n")
for (m in key_mafs) {
  sub <- summary_tab[MAF == m]
  cat(sprintf("    MAF=%.2f:\n", m))
  for (j in seq_len(nrow(sub))) {
    cat(sprintf("      %s (N=%s): min |beta| = %.3f (R² = %.1e)\n",
                sub$source[j], format(sub$N[j], big.mark = ","),
                sub$min_beta_80pct[j], sub$R2_80pct[j]))
  }
}

cat("\n  Key finding: At common variants (MAF=0.20), PanUKBB AFR requires\n")
afr_common <- results[source == "PanUKBB AFR" & MAF == 0.20]
eur_common <- results[source == "UKBB EUR" & MAF == 0.20]
cat(sprintf("    |beta| > %.3f (R² > %.1e) vs UKBB EUR |beta| > %.3f (R² > %.1e)\n",
            afr_common$min_beta_80pct, afr_common$R2_80pct,
            eur_common$min_beta_80pct, eur_common$R2_80pct))
cat(sprintf("    AFR/EUR power ratio: %.1fx fewer detectable loci\n",
            afr_common$min_beta_80pct / eur_common$min_beta_80pct))

# ==============================================================================
# 4. Save
# ==============================================================================
fwrite(results, file.path(RESULTS_DIR, "power_curves.csv"))
fwrite(summary_tab, file.path(RESULTS_DIR, "power_summary.csv"))
cat("\n  Saved: power_curves.csv, power_summary.csv\n")

# ==============================================================================
# 5. Figures
# ==============================================================================
cat("\n--- Generating figures ---\n")

# Color scheme
src_colors <- c("PanUKBB AFR" = "#E66101", "PanUKBB CSA" = "#FDB863",
                "BBJ" = "#5E3C99", "UKBB EUR" = "#2C7BB6")

# (a) Minimum detectable beta vs MAF
pdf(file.path(RESULTS_DIR, "power_min_beta.pdf"), width = 9, height = 6)
p1 <- ggplot(results, aes(x = MAF, y = min_beta_80pct, color = source)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = src_colors) +
  scale_y_log10() +
  geom_hline(yintercept = c(0.05, 0.10, 0.25), linetype = "dotted", alpha = 0.3) +
  labs(title = "Minimum Detectable Effect Size by GWAS Sample Size",
       subtitle = expression(paste("80% power at ", p < 5 %*% 10^-8,
                                    "; SE(", beta, ") ", approx, " 1/", sqrt(2*N*MAF*(1-MAF)))),
       x = "Minor Allele Frequency", y = expression(paste("Min detectable |", beta, "| (log scale)")),
       color = "GWAS Source") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
print(p1)
dev.off()
cat("  Saved: power_min_beta.pdf\n")

# (b) Variance explained (R²) vs MAF
pdf(file.path(RESULTS_DIR, "power_r2.pdf"), width = 9, height = 6)
p2 <- ggplot(results, aes(x = MAF, y = R2_80pct, color = source)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = src_colors) +
  scale_y_log10(labels = scales::scientific) +
  labs(title = "Minimum Detectable Variance Explained (R²) by Sample Size",
       subtitle = "80% power at p<5e-8; larger N detects smaller effects",
       x = "Minor Allele Frequency", y = expression(paste("Min R² (log scale)")),
       color = "GWAS Source") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
print(p2)
dev.off()
cat("  Saved: power_r2.pdf\n")

# (c) N-fold power ratio relative to UKBB EUR
ratio_dt <- merge(
  results[, .(source, MAF, min_beta_80pct)],
  results[source == "UKBB EUR", .(MAF, eur_beta = min_beta_80pct)],
  by = "MAF"
)
ratio_dt[, beta_ratio := min_beta_80pct / eur_beta]

pdf(file.path(RESULTS_DIR, "power_ratio_vs_eur.pdf"), width = 9, height = 5)
p3 <- ggplot(ratio_dt[source != "UKBB EUR"],
             aes(x = MAF, y = beta_ratio, color = source)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = src_colors) +
  labs(title = "Effect Size Ratio Relative to UKBB EUR (N=361K)",
       subtitle = "Higher = needs larger effects; PanUKBB AFR needs ~7x larger effects",
       x = "Minor Allele Frequency",
       y = expression(paste("Min |", beta, "| ratio vs UKBB EUR")),
       color = "GWAS Source") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
print(p3)
dev.off()
cat("  Saved: power_ratio_vs_eur.pdf\n")

cat("\n=== Script 50d: Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
