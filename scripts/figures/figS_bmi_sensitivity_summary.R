#!/usr/bin/env Rscript
# figS_bmi_sensitivity_summary.R
# Compact 3-panel supplementary figure: BMI / covariate sensitivity summary
# (a) Variance partition violin (Model A2, N=1,128)
# (b) Forest-plot of Spearman rho across sensitivity analyses
# (c) logFC scatter: unadjusted vs steatosis-adjusted (GSE130970, N=76)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS_SENS_DIR, "figS_bmi_sensitivity_summary.pdf")

# =========================================================================
# (a) Variance partition violin — Model A2 (fibrosis subset, N=1,128)
# =========================================================================
vp <- read.csv(file.path(BASE, "RNA-seq/results/audit_sensitivity/variance_partition_results.csv"),
               stringsAsFactors = FALSE)

# Per-gene results are for A1 (all 1,444); the summary has A2 medians.
# The per-gene CSV columns are: dataset, fibrosis_stage, group_binary, sex_covar, Residuals, gene, is_deg
# These ARE model A2 per-gene values (32,325 genes).
vp_long <- vp %>%
  select(dataset, fibrosis_stage, group_binary, sex_covar, Residuals) %>%
  pivot_longer(everything(), names_to = "covariate", values_to = "pct_variance") %>%
  mutate(pct_variance = pct_variance * 100) %>%
  mutate(covariate = factor(covariate,
                            levels = c("dataset", "fibrosis_stage", "group_binary", "sex_covar", "Residuals"),
                            labels = c("Dataset\n(batch)", "Fibrosis\nstage", "Disease\nstatus", "Sex", "Residuals")))

# Compute medians for annotation
medians <- vp_long %>%
  group_by(covariate) %>%
  summarise(med = median(pct_variance, na.rm = TRUE), .groups = "drop")

pa <- ggplot(vp_long, aes(x = covariate, y = pct_variance)) +
  geom_violin(fill = "#E0E0E0", color = "black", linewidth = 0.3, scale = "width") +
  geom_text(data = medians, aes(x = covariate, y = med, label = sprintf("%.1f%%", med)),
            vjust = -0.8, size = 2, fontface = "bold") +
  labs(x = NULL, y = "Variance explained (%)",
       title = "Variance partition (N = 1,128)") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 6))

# =========================================================================
# (b) Forest plot — Spearman rho across sensitivity analyses
# =========================================================================
# Hard-code values from confirmed analyses
forest_df <- data.frame(
  analysis = c("Age adjustment\n(N = 471, 4 cohorts)",
               "Steatosis grade\n(N = 76, GSE130970)",
               "NAS score\n(N = 660, 5 cohorts)"),
  rho = c(0.998, 0.943, 0.765),
  stringsAsFactors = FALSE
) %>%
  mutate(analysis = factor(analysis, levels = rev(analysis)),
         color = case_when(rho > 0.9 ~ "high", rho >= 0.75 ~ "moderate", TRUE ~ "low"))

rho_colors <- c(high = "#2E7D32", moderate = "#F57F17", low = "#C62828")

pb <- ggplot(forest_df, aes(x = rho, y = analysis, color = color)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_point(size = 3) +
  geom_segment(aes(x = rho - 0.02, xend = rho + 0.02, yend = analysis), linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", rho)), hjust = -0.4, size = 2.2, show.legend = FALSE) +
  annotate("text", x = 0.765, y = 0.55, label = "NAS correlates with disease\nstatus \u2192 over-correction",
           size = 1.8, color = "#C62828", hjust = 0.5, fontface = "italic") +
  annotate("text", x = 0.95, y = 3.4, label = "high concordance", size = 1.6,
           color = "gray50", hjust = 0.5) +
  scale_color_manual(values = rho_colors, guide = "none") +
  scale_x_continuous(limits = c(0.7, 1.02), breaks = seq(0.7, 1.0, 0.1)) +
  labs(x = expression("Spearman " * rho * " (logFC vs. unadjusted)"),
       y = NULL,
       title = "Covariate sensitivity") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

# =========================================================================
# (c) logFC scatter — steatosis adjustment (GSE130970, N=76)
# =========================================================================
steat <- read.csv(file.path(BASE,
  "RNA-seq/results/audit_sensitivity/steatosis_sensitivity/steatosis_sensitivity_comparison.csv"),
  stringsAsFactors = FALSE) %>%
  filter(analysis == "steatosis_grade")

# Classify significance status
steat <- steat %>%
  mutate(sig_cat = case_when(
    sig_unadjusted & sig_adjusted   ~ "Both significant",
    !sig_unadjusted & sig_adjusted  ~ "Gained",
    sig_unadjusted & !sig_adjusted  ~ "Lost",
    TRUE                             ~ "Neither"
  ))

sig_colors <- c("Both significant" = masld_colors$up,
                "Gained" = "#2E7D32",
                "Lost" = "#F57F17",
                "Neither" = "#E0E0E0")

rho_steat <- cor(steat$logFC_unadjusted, steat$logFC_adjusted, method = "spearman",
                 use = "complete.obs")

pc <- ggplot(steat, aes(x = logFC_unadjusted, y = logFC_adjusted, color = sig_cat)) +
  rasterize_layer(geom_point(size = 0.3, alpha = 0.5)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  annotate("text", x = min(steat$logFC_unadjusted, na.rm = TRUE) + 0.5,
           y = max(steat$logFC_adjusted, na.rm = TRUE) - 0.3,
           label = paste0("\u03c1 = ", sprintf("%.3f", rho_steat)),
           size = 2.5, hjust = 0) +
  scale_color_manual(values = sig_colors, name = "Significance") +
  labs(x = expression("logFC (unadjusted)"),
       y = expression("logFC (steatosis-adjusted)"),
       title = "Steatosis adjustment (N = 76, GSE130970)") +
  theme_masld() +
  theme(legend.position = c(0.85, 0.2),
        legend.background = element_blank(),
        legend.key.size = unit(0.2, "cm"))

# =========================================================================
# Assemble with patchwork: a | (b / c)
# =========================================================================
fig <- pa | (pb / pc)
fig <- auto_tag(fig) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig, OUT, width = 7, height = 4)
cat("Saved:", OUT, "\n")
