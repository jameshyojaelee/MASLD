#!/usr/bin/env Rscript
# =============================================================================
# fig_c13_fibrosis_confounding.R
# Multi-panel publication figure: C13 fibrosis-confounding discovery
# Shows that 93.5% of NASH-vs-NAFL DEGs are fibrosis-confounded
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

# --- Paths -------------------------------------------------------------------
base_dir <- Sys.getenv("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
results_dir <- file.path(base_dir, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")

# Source publication theme and centralized output paths
source(file.path(base_dir,
  "RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions/theme_publication.R"))
source(file.path(base_dir, "scripts/figures/load_figure_data.R"))

out_dir   <- FIGS02_DIR
panel_dir <- file.path(out_dir, "panels")
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

# --- Colors ------------------------------------------------------------------
col_c2       <- sanjana_colors[["Grey"]]       # #b0b0b0
col_c13      <- sanjana_colors[["Magenta"]]     # #e14b9d
col_robust   <- sanjana_colors[["Magenta"]]     # #e14b9d  (robust NASH core)
col_confused <- sanjana_colors[["Grey"]]         # #b0b0b0  (fibrosis-confounded)
col_masked   <- sanjana_colors[["Blue"]]         # #4baeef  (fibrosis-masked)

# --- Load data ---------------------------------------------------------------
cat("Loading data...\n")

# C13: fibrosis-adjusted NASH vs NAFL (padj column)
c13 <- read.csv(file.path(results_dir, "progression/c13_nash_vs_nafl_fib_adj_dream.csv"),
                 stringsAsFactors = FALSE)

# C2: unadjusted NASH vs NAFL (adj.P.Val column)
c2 <- read.csv(file.path(results_dir, "disease_signatures/nafl_vs_nash_dream.csv"),
                stringsAsFactors = FALSE)

# Harmonize p-value column names
if (!"padj" %in% names(c2) && "adj.P.Val" %in% names(c2)) {
  c2$padj <- c2$adj.P.Val
}
if (!"padj" %in% names(c13) && "adj.P.Val" %in% names(c13)) {
  c13$padj <- c13$adj.P.Val
}

# GSEA results
gsea <- read.csv(file.path(results_dir, "progression/progression_gsea_all.csv"),
                  stringsAsFactors = FALSE)

# Contrast summary
contrast_summary <- read.csv(file.path(results_dir, "progression/progression_contrast_summary.csv"),
                              stringsAsFactors = FALSE)

cat(sprintf("  C2: %d genes, %d DEGs (padj<0.1)\n", nrow(c2), sum(c2$padj < 0.1, na.rm = TRUE)))
cat(sprintf("  C13: %d genes, %d DEGs (padj<0.1)\n", nrow(c13), sum(c13$padj < 0.1, na.rm = TRUE)))

# --- Merge C2 and C13 --------------------------------------------------------
# Use gene column (Ensembl ID) to merge; symbol from c13 (both have it)
merged <- inner_join(
  c2 %>% select(gene, logFC_c2 = logFC, padj_c2 = padj),
  c13 %>% select(gene, logFC_c13 = logFC, padj_c13 = padj, symbol),
  by = "gene"
)

# Classify genes
merged <- merged %>%
  mutate(
    sig_c2  = padj_c2  < 0.1,
    sig_c13 = padj_c13 < 0.1,
    category = case_when(
      sig_c2 & sig_c13  ~ "Robust NASH core",
      sig_c2 & !sig_c13 ~ "Fibrosis-confounded",
      !sig_c2 & sig_c13 ~ "Fibrosis-masked",
      TRUE              ~ "NS in both"
    )
  )

n_robust    <- sum(merged$category == "Robust NASH core")
n_confounded <- sum(merged$category == "Fibrosis-confounded")
n_masked    <- sum(merged$category == "Fibrosis-masked")

cat(sprintf("  Robust NASH core: %d\n", n_robust))
cat(sprintf("  Fibrosis-confounded (C2-only): %d\n", n_confounded))
cat(sprintf("  Fibrosis-masked (C13-only): %d\n", n_masked))

# Compute Spearman rho
rho <- cor(merged$logFC_c2, merged$logFC_c13, method = "spearman", use = "complete.obs")
cat(sprintf("  Spearman rho (C2 vs C13 logFC): %.3f\n", rho))

# =============================================================================
# PANEL A: DEG count bar chart + breakdown
# =============================================================================

# Left side: total DEG bars
deg_totals <- data.frame(
  contrast = c("C2\n(Unadjusted)", "C13\n(Fibrosis-adjusted)"),
  count    = c(sum(merged$sig_c2), sum(merged$sig_c13 | (merged$sig_c2 & merged$sig_c13))),
  stringsAsFactors = FALSE
)
# Actually use the original DEG counts
n_c2_deg  <- sum(c2$padj < 0.1, na.rm = TRUE)
n_c13_deg <- sum(c13$padj < 0.1, na.rm = TRUE)
pct_reduction <- round((1 - n_c13_deg / n_c2_deg) * 100, 1)

# Stacked breakdown bar
breakdown_df <- data.frame(
  category = factor(
    c("Robust NASH core\n(both significant)",
      "Fibrosis-confounded\n(C2-only)",
      "Fibrosis-masked\n(C13-only)"),
    levels = c("Fibrosis-confounded\n(C2-only)",
               "Robust NASH core\n(both significant)",
               "Fibrosis-masked\n(C13-only)")
  ),
  count = c(n_robust, n_confounded, n_masked),
  stringsAsFactors = FALSE
)

# Combined panel A: side-by-side total bars + breakdown bars
bar_df <- data.frame(
  label = factor(
    c("C2\n(Unadjusted)", "C13\n(Fibrosis-adj.)"),
    levels = c("C2\n(Unadjusted)", "C13\n(Fibrosis-adj.)")
  ),
  count = c(n_c2_deg, n_c13_deg),
  fill  = c("C2", "C13"),
  stringsAsFactors = FALSE
)

pA_bars <- ggplot(bar_df, aes(x = label, y = count, fill = fill)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = formatC(count, format = "d", big.mark = ",")),
            vjust = -0.5, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("C2" = col_c2, "C13" = col_c13)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = comma) +
  # Arrow annotation for reduction
  annotate("segment",
           x = 1.05, xend = 1.95,
           y = n_c2_deg * 0.85, yend = n_c13_deg + 300,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           color = col_c13, linewidth = 0.5) +
  annotate("text", x = 1.5, y = n_c2_deg * 0.65,
           label = paste0(pct_reduction, "%\nreduction"),
           color = col_c13, size = 3, fontface = "bold") +
  labs(x = NULL, y = "DEGs (padj < 0.1)", title = "NASH vs NAFL DEGs") +
  theme_publication() +
  theme(panel.grid.major.x = element_blank())

# Breakdown horizontal bar
breakdown_df2 <- data.frame(
  category = factor(
    c("Robust NASH core", "Fibrosis-confounded", "Fibrosis-masked"),
    levels = c("Fibrosis-masked", "Robust NASH core", "Fibrosis-confounded")
  ),
  count = c(n_robust, n_confounded, n_masked),
  stringsAsFactors = FALSE
)

pA_breakdown <- ggplot(breakdown_df2, aes(x = count, y = "Breakdown", fill = category)) +
  geom_col(position = "stack", width = 0.5) +
  geom_text(aes(label = paste0(category, "\n(n=", formatC(count, big.mark = ","), ")")),
            position = position_stack(vjust = 0.5),
            size = 1.8, color = "white", fontface = "bold", lineheight = 0.85) +
  scale_fill_manual(values = c(
    "Robust NASH core"    = col_robust,
    "Fibrosis-confounded" = col_confused,
    "Fibrosis-masked"     = col_masked
  )) +
  scale_x_continuous(labels = comma) +
  labs(x = "Number of genes", y = NULL, title = "Gene classification") +
  theme_publication() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  )

pA <- pA_bars / pA_breakdown + plot_layout(heights = c(3, 1))

# =============================================================================
# PANEL B: Scatter plot C2 logFC vs C13 logFC
# =============================================================================

# Set factor levels for plotting order (NS first, then confounded, then robust, masked on top)
merged$category <- factor(merged$category,
  levels = c("NS in both", "Fibrosis-confounded", "Robust NASH core", "Fibrosis-masked"))

# Genes to label (robust core + key fibrosis-confounded)
label_genes <- c("TREM2", "LPL", "FABP4", "MMP9", "IL1B", "COL1A1",
                 "TGFB1", "FABP5", "CCL20", "IGFBP2", "SHBG",
                 "COL3A1", "VCAN", "THBS2", "TIMP1", "COL1A2",
                 "AKR1B10", "OLR1", "PBK", "ADAMDEC1")

merged_label <- merged %>%
  filter(symbol %in% label_genes & category != "NS in both")

pB <- ggplot(merged %>% filter(category != "NS in both" | runif(n()) < 0.05),
             aes(x = logFC_c2, y = logFC_c13)) +
  # NS points (thinned for speed)
  geom_point(data = . %>% filter(category == "NS in both"),
             color = "grey88", size = 0.2, alpha = 0.3) +
  # Fibrosis-confounded

  geom_point(data = . %>% filter(category == "Fibrosis-confounded"),
             color = col_confused, size = 0.6, alpha = 0.5) +
  # Robust NASH core
  geom_point(data = . %>% filter(category == "Robust NASH core"),
             color = col_robust, size = 0.8, alpha = 0.7) +
  # Fibrosis-masked
  geom_point(data = . %>% filter(category == "Fibrosis-masked"),
             color = col_masked, size = 0.8, alpha = 0.7) +
  # Diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  # Labels
  geom_text_repel(
    data = merged_label,
    aes(label = symbol),
    size = 2, fontface = "italic",
    max.overlaps = 25,
    segment.size = 0.2, segment.color = "grey50",
    min.segment.length = 0.1,
    box.padding = 0.3, point.padding = 0.15,
    color = ifelse(merged_label$category == "Robust NASH core", col_robust,
            ifelse(merged_label$category == "Fibrosis-confounded", "grey40", col_masked))
  ) +
  # Annotation
  annotate("text", x = -0.8, y = 1.7,
           label = paste0("rho = ", sprintf("%.3f", rho)),
           size = 3, hjust = 0, color = "grey30") +
  labs(
    x = expression("C2 log"[2]*"FC (Unadjusted NASH vs NAFL)"),
    y = expression("C13 log"[2]*"FC (Fibrosis-adjusted)"),
    title = expression("Effect size: unadjusted vs fibrosis-adjusted")
  ) +
  coord_cartesian(xlim = c(-2, 2.2), ylim = c(-2, 2.2)) +
  theme_publication() +
  theme(legend.position = "none")

# =============================================================================
# PANEL C: Pathway heatmap — C2 vs C13 NES comparison
# =============================================================================

# Extract Hallmark GSEA for C2 and C13
gsea_c2 <- gsea %>%
  filter(contrast_id == "C2_NASH_vs_NAFL", collection == "Hallmark") %>%
  select(pathway, NES_c2 = NES, padj_c2 = padj)

gsea_c13 <- gsea %>%
  filter(contrast_id == "C13_NASH_vs_NAFL_FibAdj", collection == "Hallmark") %>%
  select(pathway, NES_c13 = NES, padj_c13 = padj)

gsea_merged <- inner_join(gsea_c2, gsea_c13, by = "pathway")

# Clean pathway names
gsea_merged$pathway_clean <- gsub("^HALLMARK_", "", gsea_merged$pathway)
gsea_merged$pathway_clean <- gsub("_", " ", gsea_merged$pathway_clean)
gsea_merged$pathway_clean <- tools::toTitleCase(tolower(gsea_merged$pathway_clean))

# Select top 15 pathways by absolute C2 NES
top_pathways <- gsea_merged %>%
  arrange(desc(abs(NES_c2))) %>%
  slice_head(n = 15)

# Pivot for heatmap
heat_long <- top_pathways %>%
  select(pathway_clean, NES_c2, NES_c13, padj_c2, padj_c13) %>%
  pivot_longer(
    cols = c(NES_c2, NES_c13),
    names_to = "contrast",
    values_to = "NES"
  ) %>%
  mutate(
    padj = ifelse(contrast == "NES_c2", padj_c2, padj_c13),
    sig  = padj < 0.05,
    contrast_label = ifelse(contrast == "NES_c2", "C2\n(Unadjusted)", "C13\n(Fibrosis-adj.)"),
    # Mark significance
    sig_label = ifelse(sig, "", "\u2020")
  )

# Order pathways by C2 NES (descending)
pathway_order <- top_pathways %>% arrange(NES_c2) %>% pull(pathway_clean)
heat_long$pathway_clean <- factor(heat_long$pathway_clean, levels = pathway_order)
heat_long$text_color <- ifelse(abs(heat_long$NES) > 2.3, "white", "black")

pC <- ggplot(heat_long, aes(x = contrast_label, y = pathway_clean, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", NES), color = text_color),
            size = 2, show.legend = FALSE) +
  scale_color_identity() +
  # NS marker
  geom_point(data = heat_long %>% filter(!sig),
             aes(x = contrast_label, y = pathway_clean),
             shape = 4, size = 2, color = "grey30", stroke = 0.5,
             inherit.aes = FALSE) +
  scale_fill_gradient2(
    low = sanjana_colors[["Blue"]], mid = "white", high = sanjana_colors[["Magenta"]],
    midpoint = 0, limits = c(-2.5, 3.5),
    name = "NES"
  ) +
  labs(x = NULL, y = NULL,
       title = "Hallmark pathway enrichment") +
  theme_publication() +
  theme(
    axis.text.y = element_text(size = 5.5),
    panel.grid  = element_blank(),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.25, "cm")
  )

# =============================================================================
# PANEL D: Top 20 robust core genes (lollipop by |logFC|)
# =============================================================================

robust_genes <- merged %>%
  filter(category == "Robust NASH core") %>%
  arrange(desc(abs(logFC_c13))) %>%
  slice_head(n = 20) %>%
  mutate(
    direction = ifelse(logFC_c13 > 0, "Up", "Down"),
    symbol = ifelse(symbol == "" | is.na(symbol), gene, symbol)
  )

robust_genes$symbol <- factor(robust_genes$symbol,
                               levels = robust_genes$symbol[order(abs(robust_genes$logFC_c13))])

pD <- ggplot(robust_genes, aes(x = logFC_c13, y = symbol, color = direction)) +
  geom_segment(aes(x = 0, xend = logFC_c13, y = symbol, yend = symbol),
               linewidth = 0.4) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_color_manual(values = c("Up" = col_c13, "Down" = col_masked),
                     name = "Direction") +
  labs(
    x = expression("log"[2]*"FC (C13, fibrosis-adjusted)"),
    y = NULL,
    title = "Top 20 robust NASH core genes"
  ) +
  theme_publication() +
  theme(
    axis.text.y  = element_text(face = "italic", size = 6),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    panel.grid.major.y = element_blank()
  )

# =============================================================================
# Assemble composite figure
# =============================================================================

composite <- (pA | pB) / (pC | pD) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 10, face = "bold", family = font_family)
    )
  ) +
  plot_layout(heights = c(1, 1))

# --- Save outputs ------------------------------------------------------------
cat("Saving figures...\n")

# Composite
save_pdf(composite,
         file.path(out_dir, "fig_c13_fibrosis_confounding.pdf"),
         width = 12, height = 8)
cat(sprintf("  Saved: %s\n", file.path(out_dir, "fig_c13_fibrosis_confounding.pdf")))

# Individual panels → panels/ subdirectory
save_pdf(pA,
         file.path(panel_dir, "panel_c13_a_deg_reduction.pdf"),
         width = 4, height = 4)

save_pdf(pB,
         file.path(panel_dir, "panel_c13_b_scatter.pdf"),
         width = 5, height = 5)

save_pdf(pC,
         file.path(panel_dir, "panel_c13_c_pathway_heatmap.pdf"),
         width = 4.5, height = 4.5)

save_pdf(pD,
         file.path(panel_dir, "panel_c13_d_robust_core_lollipop.pdf"),
         width = 4, height = 4)

cat("Done.\n")
cat(sprintf("\nSummary:\n"))
cat(sprintf("  C2 (unadjusted):     %s DEGs\n", formatC(n_c2_deg, big.mark = ",")))
cat(sprintf("  C13 (fibrosis-adj.): %s DEGs\n", formatC(n_c13_deg, big.mark = ",")))
cat(sprintf("  Reduction:           %s%%\n", pct_reduction))
cat(sprintf("  Robust NASH core:    %d genes\n", n_robust))
cat(sprintf("  Fibrosis-confounded: %s genes\n", formatC(n_confounded, big.mark = ",")))
cat(sprintf("  Fibrosis-masked:     %d genes\n", n_masked))
cat(sprintf("  Spearman rho:        %.3f\n", rho))
