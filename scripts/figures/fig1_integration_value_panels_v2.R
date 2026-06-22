#!/usr/bin/env Rscript
# =============================================================================
# Fig 1 candidate panels — "Why Integration Matters" (batch 2)
# Generates 4 standalone PDF plots:
#   5. Per-study detection heatmap (gene × cohort significance grid)
#   6. Variance partitioning (batch vs disease vs sex vs residual)
#   7. Effect size precision gain (per-study SE vs integrated SE)
#   8. Heterogeneity vs integration power (I² vs significance)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ─── Shared constants ────────────────────────────────────────────────────────
STUDY_NAMES <- c(
  GSE213621 = "GSE213621", GSE135251 = "GSE135251", GSE130970 = "GSE130970",
  GSE162694 = "GSE162694", GSE174478 = "GSE174478", GSE193066 = "GSE193066",
  GSE240729 = "GSE240729", GSE126848 = "GSE126848"
)
COMPARE_STUDIES <- names(STUDY_NAMES)
N_STUDIES <- length(COMPARE_STUDIES)

PANEL_DIR <- file.path(FIG1_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ─── Load data ───────────────────────────────────────────────────────────────
cat("Loading data...\n")
dream     <- load_dream_results()
per_study <- load_per_study_de()
meta      <- load_meta_results()
vp        <- load_variance_partition()

dream[, gene_clean := sub("\\..*", "", gene)]
per_study[, gene_clean := sub("\\..*", "", gene)]

ps <- per_study[dataset %in% COMPARE_STUDIES]


# =============================================================================
# IDEA 5: Per-Study Detection Heatmap
# =============================================================================
cat("\n── Idea 5: Per-study detection heatmap ──\n")

# Select top integrated DEGs: 15 upregulated + 15 downregulated by |t|
sig_degs <- dream[bulk_padj < 0.01]
top_up  <- sig_degs[bulk_logFC > 0][order(-abs(t))][1:15]
top_dn  <- sig_degs[bulk_logFC < 0][order(-abs(t))][1:15]
hm_genes <- rbind(top_up, top_dn)
hm_genes <- hm_genes[!is.na(symbol) & symbol != ""]
# Remove duplicates and unknowns
hm_genes <- hm_genes[!duplicated(symbol)]
hm_genes <- hm_genes[!grepl("^ENSG", symbol)]
# Supplement if needed
if (nrow(hm_genes) < 25) {
  extra <- sig_degs[!symbol %in% hm_genes$symbol & !grepl("^ENSG", symbol)][
    order(-abs(t))][1:(30 - nrow(hm_genes))]
  hm_genes <- rbind(hm_genes, extra)
}
hm_genes <- hm_genes[seq_len(min(30, nrow(hm_genes)))]

cat(sprintf("  Selected %d genes for heatmap\n", nrow(hm_genes)))

# Build gene × study significance matrix
hm_ps <- ps[gene_clean %in% hm_genes$gene_clean,
            .(gene_clean, dataset, padj, logFC)]
hm_ps <- merge(hm_ps, hm_genes[, .(gene_clean, symbol, bulk_logFC)],
               by = "gene_clean")

# Significance status: significant + concordant direction
hm_ps[, sig_status := fifelse(
  padj < 0.1 & sign(logFC) == sign(bulk_logFC), "Significant",
  fifelse(padj < 0.1, "Opposite direction", "Not significant")
)]
hm_ps[, author := STUDY_NAMES[dataset]]

# Gene order: sorted by dream logFC (up top, down bottom)
gene_order <- hm_genes[order(-bulk_logFC), symbol]
hm_ps[, symbol := factor(symbol, levels = gene_order)]

# Study order: by sample size (largest left)
study_order <- c("GSE213621", "GSE135251", "GSE162694", "GSE240729",
                 "GSE193066", "GSE174478", "GSE130970", "GSE126848")
hm_ps[, author := factor(author, levels = study_order)]

# Add integrated result as a column
hm_dream <- hm_genes[, .(symbol, bulk_padj,
                          sig_status = fifelse(bulk_padj < 0.1,
                                               "Significant", "Not significant"))]
hm_dream[, author := factor("Integrated", levels = c(study_order, "Integrated"))]
hm_dream[, symbol := factor(symbol, levels = gene_order)]
hm_ps[, author := factor(author, levels = c(study_order, "Integrated"))]
hm_all <- rbind(
  hm_ps[, .(symbol, author, sig_status)],
  hm_dream[, .(symbol, author, sig_status)],
  fill = TRUE
)

# Fill in missing gene-study pairs as "Not tested"
full_grid <- CJ(symbol = factor(gene_order, levels = gene_order),
                author = factor(c(study_order, "Integrated"),
                                levels = c(study_order, "Integrated")))
hm_all <- merge(full_grid, hm_all, by = c("symbol", "author"), all.x = TRUE)
hm_all[is.na(sig_status), sig_status := "Not tested"]

# Colors
sig_cols <- c("Significant" = masld_colors$up,
              "Opposite direction" = "#42A5F5",
              "Not significant" = "#E0E0E0",
              "Not tested" = "white")

# Count per-gene: how many studies significant
n_per_gene <- hm_all[author != "Integrated" & sig_status == "Significant",
                      .N, by = symbol]
n_per_gene <- merge(data.table(symbol = factor(gene_order, levels = gene_order)),
                    n_per_gene, by = "symbol", all.x = TRUE)
n_per_gene[is.na(N), N := 0]

p5 <- ggplot(hm_all, aes(x = author, y = symbol, fill = sig_status)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = sig_cols, name = "Detection") +
  # Separator line before "Integrated" column
  geom_vline(xintercept = N_STUDIES + 0.5, linewidth = 0.8, color = "gray30") +
  # Gene count annotation on right margin
  geom_text(data = n_per_gene,
            aes(x = N_STUDIES + 2.2, y = symbol, label = paste0(N, "/", N_STUDIES),
                fill = NULL),
            size = 1.6, hjust = 0.5, color = "gray40") +
  annotate("text", x = N_STUDIES + 2.2, y = length(gene_order) + 0.8,
           label = "Per-study\ndetection", size = 1.6, color = "gray40",
           lineheight = 0.8) +
  coord_cartesian(clip = "off", xlim = c(0.5, N_STUDIES + 1.5)) +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        axis.text.y = element_text(face = "italic", size = 4.5),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(3, 25, 3, 3))

# fig1_idea5_detection_heatmap.pdf removed — not used in publication


# =============================================================================
# IDEA 6: Variance Partitioning (batch vs disease signal)
# =============================================================================
cat("\n── Idea 6: Variance partitioning ──\n")

if (!is.null(vp) && nrow(vp) > 0) {
  # Columns: gene, dataset, group_binary, sex, Residuals (+ symbol from loader)
  # C2 variance_partition.csv names the disease term 'group_binary' (was 'condition' pre-C2)
  vp_cols <- c("group_binary", "dataset", "sex", "Residuals")
  vp_num <- copy(vp[, c("gene", vp_cols), with = FALSE])
  for (col in vp_cols) vp_num[, (col) := as.numeric(get(col))]
  vp_long <- melt(vp_num, id.vars = "gene",
                  measure.vars = vp_cols,
                  variable.name = "source", value.name = "variance")
  # Clean names
  vp_long[, source := factor(source,
    levels = c("dataset", "group_binary", "sex", "Residuals"),
    labels = c("Cohort (batch)", "Disease status", "Sex", "Residual")
  )]
  vp_long <- vp_long[!is.na(source)]

  # Summary: median per source
  vp_summary <- vp_long[, .(median_var = median(variance, na.rm = TRUE),
                              q25 = quantile(variance, 0.25, na.rm = TRUE),
                              q75 = quantile(variance, 0.75, na.rm = TRUE)),
                          by = source]
  cat("  Variance partition medians:\n")
  for (i in seq_len(nrow(vp_summary)))
    cat(sprintf("    %s: %.1f%%\n", vp_summary$source[i],
                vp_summary$median_var[i] * 100))

  # Split into DEG vs non-DEG
  dream_sig_genes <- dream[bulk_padj < 0.1, sub("\\..*", "", gene)]
  vp_long[, gene_clean := sub("\\..*", "", gene)]
  vp_long[, is_deg := gene_clean %in% dream_sig_genes]
  vp_long[, deg_label := fifelse(is_deg,
    "Integrated DEGs\n(padj < 0.1)", "Non-significant\ngenes")]

  # Source colors
  source_cols <- c(
    "Cohort (batch)"  = "#FF8F00",    # amber
    "Disease status"   = masld_colors$up,
    "Sex"              = masld_colors$female,
    "Residual"         = "#BDBDBD"
  )

  # Box plot: variance by source, split by DEG status
  p6 <- ggplot(vp_long[source != "Residual"],
               aes(x = source, y = variance * 100, fill = source)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3,
                 linewidth = 0.3, width = 0.7) +
    facet_wrap(~ deg_label) +
    scale_fill_manual(values = source_cols, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL,
         y = "Variance explained per gene") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          strip.text = element_text(size = 6.5))

  save_fig(p6, file.path(PANEL_DIR, "fig1_idea6_variance_partition.pdf"),
           width = fig_half_width, height = 3.5)
  cat("  Saved fig1_idea6_variance_partition.pdf\n")
} else {
  cat("  WARNING: Variance partition data not available, skipping\n")
}


# =============================================================================
# IDEA 7: Effect Size Precision Gain (SE reduction)
# =============================================================================
cat("\n── Idea 7: Effect size precision gain ──\n")

# Integrated SE: |logFC / t|
dream_se <- dream[abs(t) > 0 & !is.na(bulk_padj),
                   .(gene_clean, dream_se = abs(bulk_logFC / t),
                     bulk_padj, bulk_logFC, symbol)]

# Per-study SE: compute per gene, take median across studies
ps_se <- ps[abs(t) > 0, .(gene_clean, dataset, ps_se = abs(logFC / t))]
ps_se_med <- ps_se[, .(median_ps_se = median(ps_se),
                         n_studies = .N), by = gene_clean]

# Merge
se_cmp <- merge(dream_se, ps_se_med, by = "gene_clean")
se_cmp <- se_cmp[n_studies >= 3]  # genes in ≥3 studies
se_cmp[, sig := bulk_padj < 0.1]
se_cmp[, fold_reduction := median_ps_se / dream_se]

median_fold <- median(se_cmp$fold_reduction, na.rm = TRUE)
cat(sprintf("  Genes compared: %s\n", comma(nrow(se_cmp))))
cat(sprintf("  Median SE fold-reduction: %.1fx\n", median_fold))
cat(sprintf("  Integrated DEGs: %s; Non-sig: %s\n",
            comma(sum(se_cmp$sig)), comma(sum(!se_cmp$sig))))

# Subsample for plotting (too many points otherwise)
set.seed(42)
se_plot <- se_cmp[sample(.N, min(.N, 8000))]

# Axis limits (clip outliers)
ax_max <- quantile(c(se_plot$median_ps_se, se_plot$dream_se), 0.995, na.rm = TRUE)

p7 <- ggplot(se_plot, aes(x = median_ps_se, y = dream_se)) +
  # Diagonal reference (no improvement)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.3) +
  # Points
  rasterize_layer(
    geom_point(aes(color = sig), size = 0.3, alpha = 0.4, shape = 16)
  ) +
  scale_color_manual(values = c("TRUE" = masld_colors$up,
                                 "FALSE" = "#BDBDBD"),
                     labels = c("TRUE" = "Integrated DEG",
                                "FALSE" = "Not significant"),
                     name = NULL) +
  # Fold-reduction annotation
  annotate("text", x = ax_max * 0.15, y = ax_max * 0.85,
           label = sprintf("Median SE reduction: %.1f\u00d7", median_fold),
           size = 2.5, fontface = "bold", color = masld_colors$up) +
  annotate("text", x = ax_max * 0.15, y = ax_max * 0.77,
           label = "Points below diagonal =\nmore precise integrated estimate",
           size = 1.8, color = "gray40", lineheight = 0.85) +
  coord_cartesian(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
  labs(x = "Per-study SE (median across cohorts)",
       y = "Integrated analysis SE") +
  theme_masld() +
  theme(legend.position = c(0.75, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.8),
                                          color = NA),
        legend.key.size = unit(0.25, "cm"))

save_fig(p7, file.path(PANEL_DIR, "fig1_idea7_precision_gain.pdf"),
         width = fig_half_width, height = 3.5)
cat("  Saved fig1_idea7_precision_gain.pdf\n")


# =============================================================================
# IDEA 8: I² Heterogeneity vs Integration Power
# =============================================================================
cat("\n── Idea 8: Heterogeneity vs integration power ──\n")

if (!is.null(meta) && nrow(meta) > 0) {
  # Merge meta I² with dream results
  meta[, gene_clean := sub("\\..*", "", gene)]
  het <- merge(dream[, .(gene_clean, bulk_padj, bulk_logFC, symbol)],
               meta[, .(gene_clean, meta_I2, meta_padj)],
               by = "gene_clean")
  het <- het[!is.na(meta_I2) & !is.na(bulk_padj)]

  # Classification
  het[, status := fifelse(
    bulk_padj < 0.1 & meta_padj < 0.1, "Both significant",
    fifelse(bulk_padj < 0.1, "Integrated only",
    fifelse(meta_padj < 0.1, "Meta only", "Neither"))
  )]

  cat(sprintf("  Genes with I² data: %s\n", comma(nrow(het))))
  cat("  Status counts:\n")
  print(het[, .N, by = status][order(-N)])

  # Subsample for plotting
  set.seed(42)
  het_plot <- het[sample(.N, min(.N, 10000))]

  # Known genes to label (high I², significant in dream)
  label_genes <- c("TREM2", "SPP1", "COL1A1", "CYP7A1", "CIDEC", "FAP",
                   "ACTA2", "GDF15", "TIMP1", "THY1", "PLIN2", "SLC27A5")
  het_labels <- het[symbol %in% label_genes & bulk_padj < 0.1]
  het_labels <- het_labels[!duplicated(symbol)]

  status_cols <- c(
    "Both significant"  = masld_colors$up,
    "Integrated only"   = "#7B1FA2",
    "Meta only"         = "#42A5F5",
    "Neither"           = "#E0E0E0"
  )

  p8 <- ggplot(het_plot, aes(x = meta_I2, y = -log10(bulk_padj))) +
    # Background quadrants
    annotate("rect", xmin = 50, xmax = 100, ymin = -log10(0.1), ymax = Inf,
             fill = "#FCE4EC", alpha = 0.3) +
    annotate("text", x = 75, y = max(-log10(het_plot$bulk_padj)) * 0.95,
             label = "High heterogeneity\nresolved by integration",
             size = 1.8, color = "#880E4F", fontface = "italic",
             lineheight = 0.85) +
    # Significance threshold
    geom_hline(yintercept = -log10(0.1), linetype = "dashed",
               color = "gray50", linewidth = 0.3) +
    annotate("text", x = 2, y = -log10(0.1),
             label = "padj = 0.1", hjust = 0, vjust = -0.3,
             size = 1.8, color = "gray50") +
    # Points
    rasterize_layer(
      geom_point(aes(color = status), size = 0.4, alpha = 0.4, shape = 16)
    ) +
    scale_color_manual(values = status_cols, name = "Significance") +
    # Label known genes
    geom_point(data = het_labels,
               aes(x = meta_I2, y = -log10(bulk_padj)),
               shape = 1, size = 2, color = "#880E4F", stroke = 0.4) +
    geom_text_repel(data = het_labels,
                    aes(x = meta_I2, y = -log10(bulk_padj), label = symbol),
                    size = 1.8, fontface = "italic", color = "#880E4F",
                    segment.size = 0.2, box.padding = 0.2,
                    max.overlaps = 20, seed = 42) +
    scale_x_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 100)) +
    labs(x = expression(I^2~"(inter-study heterogeneity)"),
         y = expression(-log[10]~"(integrated padj)")) +
    theme_masld() +
    theme(legend.position = c(0.18, 0.82),
          legend.background = element_rect(fill = alpha("white", 0.85),
                                            color = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 5))

  save_fig(p8, file.path(PANEL_DIR, "fig1_idea8_heterogeneity_power.pdf"),
           width = fig_half_width, height = 3.8)
  cat("  Saved fig1_idea8_heterogeneity_power.pdf\n")
} else {
  cat("  WARNING: Meta-analysis results not available, skipping\n")
}


cat("\n=== All panels generated in figures/panels/ ===\n")
cat("  5. (fig1_idea5_detection_heatmap removed)\n")
cat("  6. fig1_idea6_variance_partition.pdf\n")
cat("  7. fig1_idea7_precision_gain.pdf\n")
cat("  8. fig1_idea8_heterogeneity_power.pdf\n")
