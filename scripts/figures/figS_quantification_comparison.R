#!/usr/bin/env Rscript
# figS_quantification_comparison.R — Supplementary figure: STAR vs Kallisto
# quantification comparison supporting the switch from STAR+featureCounts
# to Kallisto pseudoalignment for the human dream mega-analysis.
#
# KEY MESSAGE: Kallisto reduces cross-cohort heterogeneity (I² 90.3→65.4%,
# τ² 76% reduction) while preserving effect-size concordance (r=0.894),
# detecting 14% more genes, and producing a more conservative DEG set.
#
# Output: figures/supplementary/figS_methods_validation/quantification/figS_quantification_comparison.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

INT_DIR  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                      "results/integration")
COMP_DIR <- file.path(BASE, "RNA-seq/results/kallisto_vs_star_comparison",
                      "kallisto_vs_star_comparison/tables")

# ── Load data ────────────────────────────────────────────────────────────────

merged         <- fread(file.path(COMP_DIR, "dream_lfc_padj_merged.csv"))
per_gene_rho   <- fread(file.path(COMP_DIR, "per_gene_spearman_log_tpm.csv"))

kall      <- fread(file.path(INT_DIR, "dream_results_ashr.csv"))
star      <- fread(file.path(INT_DIR, "dream_results_ashr_star.csv"))
meta_kall <- fread(file.path(INT_DIR, "meta_analysis_results.csv"))
meta_star <- fread(file.path(INT_DIR, "meta_analysis_results_moderated_SE_backup.csv"))

# ── Colours — method comparison palette ──────────────────────────────────────

col_star      <- "#9E9E9E"   # gray (legacy method, receding)
col_kall      <- masld_colors$up   # #C9265E magenta (focal method)
col_both      <- masld_colors$down # #1565C0 blue
col_star_only <- "#F4A674"         # warm peach
col_kall_only <- "#00695C"         # teal

# ── Panel A: LFC scatter ─────────────────────────────────────────────────────
# KEY MESSAGE: Effect sizes are highly concordant (r=0.894); the switch
# preserves biological signal.

PADJ_CUT <- 0.05
LFC_CUT  <- 0.5

merged[, sig_star := padj_star < PADJ_CUT & abs(logFC_star) > LFC_CUT]
merged[, sig_kall := padj_kall < PADJ_CUT & abs(logFC_kall) > LFC_CUT]
merged[, deg_class := fcase(
  sig_star & sig_kall,   "Both",
  sig_star & !sig_kall,  "STAR only",
  !sig_star & sig_kall,  "Kallisto only",
  default = "Neither"
)]
merged[, deg_class := factor(deg_class,
                             levels = c("Neither", "Both", "STAR only", "Kallisto only"))]

r_val  <- cor(merged$logFC_star, merged$logFC_kall, use = "complete.obs")
rho_val <- cor(merged$logFC_star, merged$logFC_kall, use = "complete.obs", method = "spearman")

n_both      <- sum(merged$deg_class == "Both")
n_star_only <- sum(merged$deg_class == "STAR only")
n_kall_only <- sum(merged$deg_class == "Kallisto only")

pA <- ggplot(merged, aes(logFC_star, logFC_kall)) +
  rasterize_layer(
    geom_point(data = merged[deg_class == "Neither"],
               color = "#E8E8E8", size = 0.2, alpha = 0.4, shape = 16)
  ) +
  rasterize_layer(
    geom_point(data = merged[deg_class != "Neither"],
               aes(color = deg_class), size = 0.45, alpha = 0.8, shape = 16)
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.25) +
  scale_color_manual(
    values = c("Both" = col_both, "STAR only" = col_star_only,
               "Kallisto only" = col_kall_only),
    guide = "none"
  ) +
  annotate("text", x = -3.2, y = 3.1,
           label = sprintf("r = %.3f  |  ρ = %.3f", r_val, rho_val),
           size = GEOM_TEXT_6PT, hjust = 0, vjust = 1, color = "black", fontface = "plain") +
  annotate("text", x = 2.8, y = -2.5,
           label = sprintf("Both: %s\nSTAR only: %s\nKallisto only: %s",
                           comma(n_both), comma(n_star_only), comma(n_kall_only)),
           size = GEOM_TEXT_6PT, hjust = 1, vjust = 0, color = "black") +
  coord_fixed(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5)) +
  labs(x = expression(log[2]*"FC (STAR + featureCounts)"),
       y = expression(log[2]*"FC (Kallisto + tximport)")) +
  theme_masld(base_size = 7) +
  theme(plot.margin = margin(4, 4, 2, 4))

# ── Panel B: DEG overlap ─────────────────────────────────────────────────────
# KEY MESSAGE: Kallisto is more conservative; 626 STAR-only DEGs are likely
# quantification artifacts (multimapper / length-bias).

jaccard <- n_both / (n_both + n_star_only + n_kall_only)

deg_counts <- data.table(
  category = factor(c("STAR\nonly", "Shared", "Kallisto\nonly"),
                    levels = c("STAR\nonly", "Shared", "Kallisto\nonly")),
  count = c(n_star_only, n_both, n_kall_only),
  col   = c("star_only", "both", "kall_only")
)

pB <- ggplot(deg_counts, aes(category, count, fill = col)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = comma(count)), vjust = -0.4, size = GEOM_TEXT_6PT) +
  annotate("text", x = 2, y = max(deg_counts$count) * 0.85,
           label = sprintf("Jaccard = %.3f", jaccard),
           size = GEOM_TEXT_6PT, fontface = "plain", color = "white") +
  scale_fill_manual(values = c("star_only" = col_star_only, "both" = col_both,
                               "kall_only" = col_kall_only), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Tier 1 DEGs") +
  theme_masld(base_size = 7) +
  theme(plot.margin = margin(4, 4, 2, 4),
        panel.grid.major.y = element_line(linewidth = 0.15, color = "gray92"))

# ── Panel C: I² density ──────────────────────────────────────────────────────
# KEY MESSAGE: Kallisto dramatically reduces between-cohort heterogeneity
# (median I² drops 25 percentage points).

het_kall <- meta_kall[, .(I2 = meta_I2, method = "Kallisto")]
het_star <- meta_star[, .(I2 = meta_I2, method = "STAR")]
het_all  <- rbind(het_kall, het_star)
het_all[, method := factor(method, levels = c("STAR", "Kallisto"))]

med_kall_i2 <- median(het_kall$I2, na.rm = TRUE)
med_star_i2 <- median(het_star$I2, na.rm = TRUE)

pC <- ggplot(het_all, aes(I2, fill = method, color = method)) +
  geom_density(alpha = 0.3, linewidth = 0.35) +
  geom_vline(xintercept = med_star_i2, linetype = "dashed",
             color = col_star, linewidth = 0.25) +
  geom_vline(xintercept = med_kall_i2, linetype = "dashed",
             color = col_kall, linewidth = 0.25) +
  annotate("label", x = 55, y = 0.058,
           label = sprintf("STAR: %.1f%%", med_star_i2),
           size = GEOM_TEXT_6PT, color = "black", fill = "white",
           label.padding = unit(0.12, "lines")) +
  annotate("label", x = 55, y = 0.050,
           label = sprintf("Kallisto: %.1f%%", med_kall_i2),
           size = GEOM_TEXT_6PT, color = col_kall, fill = "white",
           label.padding = unit(0.12, "lines")) +
  scale_fill_manual(values = c("STAR" = col_star, "Kallisto" = col_kall), guide = "none") +
  scale_color_manual(values = c("STAR" = "gray50", "Kallisto" = col_kall), guide = "none") +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  labs(x = expression(I^2 ~ "(%)"),
       y = "Density") +
  theme_masld(base_size = 7) +
  theme(plot.margin = margin(4, 4, 2, 4))

# ── Panel D: τ² boxplot ──────────────────────────────────────────────────────
# KEY MESSAGE: Between-cohort variance drops 76% — Kallisto's EM allocation
# and length-aware aggregation produce more consistent quantification.

tau_kall <- meta_kall[, .(tau2 = meta_tau2, method = "Kallisto")]
tau_star <- meta_star[, .(tau2 = meta_tau2, method = "STAR")]
tau_all  <- rbind(tau_kall, tau_star)
tau_all[, method := factor(method, levels = c("STAR", "Kallisto"))]

tau_meds <- tau_all[, .(med = median(tau2, na.rm = TRUE)), by = method]
pct_reduction <- 100 * (tau_meds[method == "STAR", med] - tau_meds[method == "Kallisto", med]) /
                        tau_meds[method == "STAR", med]

pD <- ggplot(tau_all, aes(method, tau2, fill = method)) +
  geom_boxplot(outlier.size = 0.15, outlier.alpha = 0.2,
               width = 0.45, linewidth = 0.25, median.linewidth = 0.6) +
  geom_text(data = tau_meds,
            aes(method, med, label = sprintf("%.4f", med)),
            vjust = -0.6, size = GEOM_TEXT_6PT, fontface = "plain") +
  annotate("text", x = 1.5, y = quantile(tau_all$tau2, 0.94, na.rm = TRUE),
           label = sprintf("%.0f%% reduction", pct_reduction),
           size = GEOM_TEXT_6PT, fontface = "plain", color = col_kall) +
  scale_fill_manual(values = c("STAR" = col_star, "Kallisto" = col_kall), guide = "none") +
  scale_y_continuous(limits = c(0, quantile(tau_all$tau2, 0.95, na.rm = TRUE)),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = expression(tau^2)) +
  theme_masld(base_size = 7) +
  theme(plot.margin = margin(4, 4, 2, 4))

# ── Panel E: Per-gene TPM concordance ─────────────────────────────────────────
# KEY MESSAGE: 64% of genes have Spearman ρ > 0.8 for per-sample TPM;
# the methods agree at the expression-level, not just the DEG level.

med_rho  <- median(per_gene_rho$spearman_log_tpm, na.rm = TRUE)
pct_high <- 100 * mean(per_gene_rho$spearman_log_tpm > 0.8, na.rm = TRUE)

pE <- ggplot(per_gene_rho, aes(spearman_log_tpm)) +
  geom_histogram(bins = 80, fill = col_kall, color = "white",
                 linewidth = 0.08, alpha = 0.85) +
  geom_vline(xintercept = med_rho, linetype = "dashed",
             color = "gray30", linewidth = 0.25) +
  annotate("label", x = 0.4, y = Inf,
           label = sprintf("median ρ = %.3f\n%.1f%% genes ρ > 0.8", med_rho, pct_high),
           size = GEOM_TEXT_6PT, hjust = 0.5, vjust = 1.3, color = "black",
           fill = "white", label.padding = unit(0.15, "lines")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     labels = label_comma()) +
  labs(x = "Per-gene Spearman ρ (log1p TPM)",
       y = "Genes") +
  theme_masld(base_size = 7) +
  theme(plot.margin = margin(4, 4, 2, 4))

# ── Panel F: Gene detection expansion ────────────────────────────────────────
# KEY MESSAGE: Kallisto detects 14% more genes via EM multimapper allocation,
# recovering lowly-expressed genes that STAR discards.

n_star_genes <- nrow(star)
n_kall_genes <- nrow(kall)
n_new        <- n_kall_genes - n_star_genes

det_dt <- data.table(
  method = factor(c("STAR", "Kallisto"), levels = c("STAR", "Kallisto")),
  count  = c(n_star_genes, n_kall_genes)
)

pF <- ggplot(det_dt, aes(method, count, fill = method)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = comma(count)), vjust = -0.4, size = GEOM_TEXT_6PT) +
  annotate("segment", x = 1.15, xend = 1.85,
           y = n_star_genes, yend = n_star_genes,
           linetype = "dotted", color = "gray50", linewidth = 0.25) +
  annotate("text", x = 2, y = (n_star_genes + n_kall_genes) / 2,
           label = sprintf("+%s\n(+%.1f%%)", comma(n_new), 100 * n_new / n_star_genes),
           size = GEOM_TEXT_6PT, fontface = "plain", color = "white") +
  scale_fill_manual(values = c("STAR" = col_star, "Kallisto" = col_kall), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = label_comma()) +
  labs(x = NULL, y = "Genes tested") +
  theme_masld(base_size = 7) +
  theme(plot.margin = margin(4, 4, 2, 4))

# ── Assemble composite ───────────────────────────────────────────────────────

fig <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 6, face = "plain"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

save_fig(fig, file.path(FIGS_QUANT_DIR, "figS_quantification_comparison.pdf"),
         width = fig_full_width, height = 7.5)
message("Saved: ", file.path(FIGS_QUANT_DIR, "figS_quantification_comparison.pdf"))

# ── Summary table ────────────────────────────────────────────────────────────

summary_dt <- data.table(
  Metric = c("Genes tested", "Tier 1 DEGs (padj<0.05, |LFC|>0.5)",
             "  Upregulated", "  Downregulated",
             "LFC Pearson r", "LFC Spearman rho",
             "DEG Jaccard", "Per-gene TPM rho (median)",
             "I2 median (%)", "tau2 median"),
  STAR = c(comma(n_star_genes),
           comma(sum(merged$sig_star)),
           comma(sum(star$padj < 0.05 & star$logFC > 0.5, na.rm = TRUE)),
           comma(sum(star$padj < 0.05 & star$logFC < -0.5, na.rm = TRUE)),
           NA_character_, NA_character_, NA_character_, NA_character_,
           sprintf("%.1f", med_star_i2),
           sprintf("%.4f", tau_meds[method == "STAR", med])),
  Kallisto = c(comma(n_kall_genes),
               comma(sum(merged$sig_kall)),
               comma(sum(kall$padj < 0.05 & kall$logFC > 0.5, na.rm = TRUE)),
               comma(sum(kall$padj < 0.05 & kall$logFC < -0.5, na.rm = TRUE)),
               NA_character_, NA_character_, NA_character_, NA_character_,
               sprintf("%.1f", med_kall_i2),
               sprintf("%.4f", tau_meds[method == "Kallisto", med])),
  Comparison = c(sprintf("+%s (+%.1f%%)", comma(n_new), 100 * n_new / n_star_genes),
                 sprintf("-%s", comma(sum(merged$sig_star) - sum(merged$sig_kall))),
                 "", "",
                 sprintf("%.3f", r_val), sprintf("%.3f", rho_val),
                 sprintf("%.3f", jaccard), sprintf("%.3f", med_rho),
                 sprintf("-%.1f pp", med_star_i2 - med_kall_i2),
                 sprintf("-%.0f%%", pct_reduction))
)
fwrite(summary_dt, file.path(FIGS_QUANT_DIR, "quantification_comparison_summary.csv"))
message("Saved: ", file.path(FIGS_QUANT_DIR, "quantification_comparison_summary.csv"))
