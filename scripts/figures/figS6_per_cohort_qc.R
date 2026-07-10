# Figure S6: Per-Cohort QC and DE
# Panels: (a) DEG counts bar, (b) Volcano grid (6 facets, rasterized)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIG_OUT, "supplementary", "archive", "figS6_per_cohort_qc.pdf")

psd <- load_per_study_de()

# T2.9 fix (2026-04-22): apply consistent padj<0.05 threshold across per-study
# DEG counts, matching `config/pipeline_params.yaml::de.padj_sig: 0.05`
# (primary significance threshold used by the dream mega-analysis and
# downstream DEG claims). The prior padj<0.1 threshold (`de.padj_de`) is
# only appropriate for exploratory annotation, not volcano-panel DEG labels.
# FIXME: the task note (T2.9) references `scripts/figures/10_volcano_plots.R:46-48`
# which does not exist in this repo; figS6_per_cohort_qc.R is the only
# per-cohort volcano figure with the same padj<0.1 bug, so the fix is applied
# here. If a future 10_volcano_plots.R is added, it should load
# PADJ_SIG from the YAML via load_figure_data.R (no new constant needed —
# just read config directly).
PADJ_SIG_CUTOFF <- 0.05  # config/pipeline_params.yaml::de.padj_sig

# ---- Panel (a): DEG counts per cohort ----
if (!is.null(psd) && nrow(psd) > 0) {
  deg_counts <- psd[padj < PADJ_SIG_CUTOFF, .N, by = dataset]
  deg_up   <- psd[padj < PADJ_SIG_CUTOFF & logFC > 0, .N, by = dataset]
  deg_down <- psd[padj < PADJ_SIG_CUTOFF & logFC < 0, .N, by = dataset]

  setnames(deg_up, "N", "Up")
  setnames(deg_down, "N", "Down")
  deg_bar <- merge(deg_up, deg_down, by = "dataset", all = TRUE)
  deg_bar[is.na(Up), Up := 0L]
  deg_bar[is.na(Down), Down := 0L]
  deg_melt <- melt(deg_bar, id.vars = "dataset", variable.name = "Direction", value.name = "Count")
  deg_melt[, Direction := factor(Direction, levels = c("Up", "Down"))]

  p_a <- ggplot(deg_melt, aes(x = reorder(dataset, -Count), y = Count, fill = Direction)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c(Up = masld_colors$up, Down = masld_colors$down), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL, y = sprintf("DEGs (padj < %s)", format(PADJ_SIG_CUTOFF))) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          legend.position = c(0.85, 0.85), legend.background = element_blank())
} else {
  p_a <- placeholder("Panel a: per-study DE data not found")
}

# ---- Panel (b): Volcano grid, 6 facets (rasterized) ----
if (!is.null(psd) && nrow(psd) > 0) {
  psd[, neg_log10p := pmin(-log10(padj), 50)]
  psd[, sig_cat := fifelse(padj < PADJ_SIG_CUTOFF & logFC >= 0.5, "Up",
                   fifelse(padj < PADJ_SIG_CUTOFF & logFC <= -0.5, "Down", "NS"))]

  # Top 3 genes per dataset for labeling
  top_genes <- psd[sig_cat != "NS", .SD[order(padj)][1:min(3, .N)], by = dataset]

  p_b <- ggplot(psd, aes(x = logFC, y = neg_log10p, color = sig_cat)) +
    rasterize_layer(geom_point(size = 0.15, alpha = 0.3, shape = 16)) +
    scale_color_manual(values = c(Up = masld_colors$up, Down = masld_colors$down,
                                  NS = masld_colors$ns), guide = "none") +
    geom_hline(yintercept = -log10(PADJ_SIG_CUTOFF),
               linetype = "dashed", linewidth = 0.2, color = "gray50") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "gray50") +
    geom_text_repel(data = top_genes, aes(label = symbol), size = GEOM_TEXT_6PT,
                    max.overlaps = 8, segment.size = 0.15, color = "black") +
    facet_wrap(~dataset, scales = "free", ncol = 3) +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~padj)) +
    theme_masld() +
    theme(strip.text = element_text(size = 6))
} else {
  p_b <- placeholder("Panel b: per-study DE data not found")
}

# ---- Assemble ----
figS6 <- p_a / p_b + plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "plain"))

save_fig_tall(figS6, OUT, height = 8)
message("[caption] Per-cohort DEG counts")
message("FigS6 saved to ", OUT)
