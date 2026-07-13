#!/usr/bin/env Rscript
# render_figS_meta_concordance.R
# Render two panels for figS_mega_validation:
#   panelD_meta_dream_scatter.pdf  — dream vs metafor logFC scatter, coloured by I²
#   panelG_meta_heterogeneity.pdf  — I² distribution + I² vs |LFC| scatter
# Input: pre-computed dream_vs_meta_concordance.csv + meta_analysis_results.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR   <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/mega_validation/panels")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))

# ── data ────────────────────────────────────────────────────────────────────
conc <- fread(file.path(RDIR, "dream_vs_meta_concordance.csv"))
meta <- fread(file.path(RDIR, "meta_analysis_results.csv"))

# summary stats for labels
n_both      <- nrow(conc)
rho_val     <- cor(conc$dream_logFC, conc$meta_logFC, method = "spearman", use = "complete.obs")  # C2-OK-sensitivity
dream_degs  <- conc[dream_sig == TRUE, .N]  # C2-OK-sensitivity
meta_degs   <- conc[meta_sig  == TRUE, .N]
both_sig    <- conc[dream_sig == TRUE & meta_sig == TRUE, .N]  # C2-OK-sensitivity
union_n     <- conc[dream_sig == TRUE | meta_sig == TRUE, .N]  # C2-OK-sensitivity
jaccard_val <- both_sig / union_n
dir_conc    <- conc[dream_sig == TRUE & meta_sig == TRUE,  # C2-OK-sensitivity
                    mean(direction_concordant, na.rm = TRUE)]

# ── Panel D: scatter logFC coloured by I² ───────────────────────────────────
lim <- max(abs(c(conc$dream_logFC, conc$meta_logFC)), na.rm = TRUE) * 1.05  # C2-OK-sensitivity

pD <- ggplot(conc, aes(dream_logFC, meta_logFC, colour = meta_I2)) +  # C2-OK-sensitivity
  geom_point(size = 0.25, alpha = 0.45) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "#B0BEC5", linewidth = 0.4) +
  scale_colour_viridis_c(
    name = expression(I^2 ~ "(%)"),
    limits = c(0, 100),
    option = "plasma", direction = -1
  ) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    x = "Dream log2FC (fixed-effect pooled)",
    y = "metafor log2FC",
    title = sprintf("Spearman rho = %.3f  |  N = %s genes",
                    rho_val, format(n_both, big.mark = ","))
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = "right",
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(0.25, "cm"))

ggsave(file.path(OUTDIR, "panelD_meta_dream_scatter.pdf"),
       pD, width = 6, height = 5.5, useDingbats = FALSE)
cat("Wrote panelD_meta_dream_scatter.pdf\n")

# ── Panel G: I² distribution + I² vs |LFC| ──────────────────────────────────
med_i2      <- median(meta$meta_I2, na.rm = TRUE)
pct_high_i2 <- 100 * mean(meta$meta_I2 > 50, na.rm = TRUE)
med_i2_degs <- median(
  meta[meta_padj < 0.05 & abs(meta_logFC) > 0.3, meta_I2], na.rm = TRUE)

pG1 <- ggplot(meta, aes(meta_I2)) +
  geom_histogram(bins = 50, fill = "#5E4FA2", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = c(25, 50, 75),
             linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  annotate("text", x = c(12.5, 37.5, 62.5, 87.5), y = Inf,
           label = c("Low", "Moderate", "Substantial", "Considerable"),
           vjust = 1.5, size = 2.2, colour = "grey40") +
  labs(
    x = expression(I^2 ~ "(%)"),
    y = "Genes",
    title = sprintf("All genes: median I² = %.1f%%  |  %.1f%% genes > 50%%",
                    med_i2, pct_high_i2),
    subtitle = sprintf("DEGs (padj<0.05 |LFC|>0.3): median I² = %.1f%%", med_i2_degs)
  ) +
  theme_masld(base_size = 7)

pG2 <- ggplot(meta, aes(abs(meta_logFC), meta_I2)) +
  geom_point(aes(colour = meta_padj < 0.05 & abs(meta_logFC) > 0.3),
             size = 0.2, alpha = 0.35) +
  scale_colour_manual(
    values = c("TRUE" = "#D6604D", "FALSE" = "#CFD8DC"),
    labels = c("TRUE" = "DEG (padj<0.05, |LFC|>0.3)", "FALSE" = "NS"),
    name = NULL
  ) +
  labs(
    x = "|Meta log2FC|",
    y = expression(I^2 ~ "(%)"),
    title = "Heterogeneity vs effect size"
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6))

pG <- pG1 + pG2 + plot_layout(ncol = 2)

ggsave(file.path(OUTDIR, "panelG_meta_heterogeneity.pdf"),
       pG, width = 10, height = 4.5, useDingbats = FALSE)
cat("Wrote panelG_meta_heterogeneity.pdf\n")

cat(sprintf(
  "\nSummary:\n  N genes in concordance table: %s\n  Spearman rho: %.3f\n",
  format(n_both, big.mark = ","), rho_val))
cat(sprintf(
  "  Dream DEGs (padj<0.1): %s | Meta DEGs: %s\n  Jaccard: %.3f | Direction concordance: %.1f%%\n",
  format(dream_degs, big.mark = ","), format(meta_degs, big.mark = ","),
  jaccard_val, 100 * dir_conc))
