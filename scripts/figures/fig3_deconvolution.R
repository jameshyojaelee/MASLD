#!/usr/bin/env Rscript
# ==========================================================================
# Figure 3: Population-scale deconvolution reveals hepatocyte-intrinsic signal
# 4-panel layout (2 rows x 2 columns):
#   a: Unadjusted vs composition-adjusted logFC scatter (RASTERIZED)
#   b: Attribution category bar chart (hep-intrinsic / comp-driven / unmasked)
#   c: BayesPrism vs MuSiC cross-method validation
#   d: Top 20 hepatocyte-intrinsic genes (lollipop)
# ==========================================================================

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

dir.create(file.path(FIGS03_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(FIGS03_DIR, "panels", "fig3_deconvolution.pdf")

# Deconvolution-specific paths
BP_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                    "results/deconvolution/bayesprism")

# ==========================================================================
# Load deconvolution data
# ==========================================================================
deconv <- load_deconv_attribution()

# Pre-initialize all panels with placeholders
p_a <- placeholder("Panel a: deconv_attribution_scores.csv not found")
p_b <- placeholder("Panel b: data unavailable")
p_c <- placeholder("Panel c: cross-method comparison not found")
p_d <- placeholder("Panel d: data unavailable")
# ==========================================================================
# Column normalization (shared across panels a, b, d)
# ==========================================================================
if (!is.null(deconv)) {
  # Detect and normalize column names
  if (!all(c("logFC_unadj", "logFC_adj", "category") %in% names(deconv))) {
    unadj_col <- intersect(c("logFC_unadj", "logFC_unadjusted", "logFC"), names(deconv))[1]
    adj_col   <- intersect(c("logFC_adj", "logFC_adjusted", "logFC_deconv"), names(deconv))[1]
    cat_col   <- intersect(c("category", "attribution_class", "class"), names(deconv))[1]
    if (!is.na(unadj_col)) setnames(deconv, unadj_col, "logFC_unadj")
    if (!is.na(adj_col))   setnames(deconv, adj_col,   "logFC_adj")
    if (!is.na(cat_col))   setnames(deconv, cat_col,   "category")
  }

  has_cols <- all(c("logFC_unadj", "logFC_adj", "category") %in% names(deconv))
} else {
  has_cols <- FALSE
}

# Full palette including Not_significant background
full_attrib_colors <- c(attribution_colors, Not_significant = masld_colors$ns)

# ==========================================================================
# Panel a: Unadjusted vs composition-adjusted logFC scatter (RASTERIZED)
# ==========================================================================
if (has_cols) {
  scatter_dt <- deconv[!is.na(logFC_unadj) & !is.na(logFC_adj)]

  # Axis limits from 99th percentile
  lim <- quantile(abs(c(scatter_dt$logFC_unadj, scatter_dt$logFC_adj)),
                  0.99, na.rm = TRUE)
  lim <- max(lim, 0.5)

  # Draw Not_significant first (background), significant categories on top
  scatter_dt[, plot_order := fifelse(category == "Not_significant", 0L, 1L)]
  setorder(scatter_dt, plot_order)

  # Dynamic legend labels with counts
  df_counts_a <- scatter_dt[, .N, by = category]
  legend_labels_a <- setNames(
    paste0(df_counts_a$category, " (n=", format(df_counts_a$N, big.mark = ","), ")"),
    df_counts_a$category
  )

  # Select top genes per significant category for labeling
  sig_scatter <- scatter_dt[category != "Not_significant" & !grepl("^ENS(G|MUSG)", symbol)]
  top_labels_a <- sig_scatter[order(-abs(logFC_adj - logFC_unadj))][1:min(12, .N)]
  if (!"symbol" %in% names(top_labels_a)) {
    top_labels_a[, symbol := gene]
  }

  p_a <- ggplot(scatter_dt,
                aes(x = logFC_unadj, y = logFC_adj, color = category)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, color = "gray50") +
    rasterize_layer(
      geom_point(size = 0.3, alpha = 0.4, shape = 16)
    ) +
    geom_point(data = top_labels_a, aes(x = logFC_unadj, y = logFC_adj),
               color = "black", size = 1, shape = 21, fill = "yellow", stroke = 0.5) +
    geom_label_repel(data = top_labels_a,
                     aes(label = symbol),
                     size = 1.6, max.overlaps = 50,
                     label.padding = 0.15, box.padding = 0.5,
                     segment.size = 0.15, show.legend = FALSE) +
    scale_color_manual(values = full_attrib_colors, name = "Attribution",
                       labels = legend_labels_a) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = expression("Unadjusted log"[2]*"FC"),
         y = expression("Composition-adjusted log"[2]*"FC")) +
    theme_masld() +
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))
}

# ==========================================================================
# Panel b: Attribution class bar chart (83% / ~8% / ~8% of sig genes)
# ==========================================================================
if (has_cols) {
  sig_deconv <- deconv[category != "Not_significant"]
  class_counts <- sig_deconv[, .N, by = category]
  total_sig <- sum(class_counts$N)
  class_counts[, pct := sprintf("%.0f%%", N / total_sig * 100)]
  class_counts[, label := gsub("_", "\n", category)]

  p_b <- ggplot(class_counts,
                aes(x = reorder(label, -N), y = N, fill = category)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(format(N, big.mark = ","), "\n(", pct, ")")),
              vjust = -0.3, size = 2.5) +
    scale_fill_manual(values = attribution_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL,
         y = paste0("Genes (N=", format(total_sig, big.mark = ","),
                    " significant)")) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 6))
}

# ==========================================================================
# Panel c: BayesPrism vs MuSiC cross-method validation
# ==========================================================================
# Try the centralized loader first
xmethod <- load_attribution_comparison()

# Fallback to the known BayesPrism directory
if (is.null(xmethod)) {
  xmethod_path <- file.path(BP_DIR, "attribution_cross_method_comparison.csv")
  if (file.exists(xmethod_path)) xmethod <- fread(xmethod_path)
}

if (!is.null(xmethod) && all(c("bp_class", "music_class") %in% names(xmethod))) {
  # Build confusion-style contingency for significant genes
  # Map BP short names to standard names
  xmethod[, bp_cat := fcase(
    bp_class == "Hepatocyte_intrinsic", "Hepatocyte_intrinsic",
    bp_class == "Composition_driven",   "Composition_driven",
    bp_class == "Unmasked",             "Unmasked",
    default = "Not_significant"
  )]
  xmethod[, music_cat := fcase(
    music_class == "Hepatocyte_intrinsic", "Hepatocyte_intrinsic",
    music_class == "Composition_driven",   "Composition_driven",
    music_class == "Unmasked",             "Unmasked",
    default = "Not_significant"
  )]

  # Compute cross-method concordance metrics
  sig_either <- xmethod[bp_cat != "Not_significant" | music_cat != "Not_significant"]
  sig_both   <- xmethod[bp_cat != "Not_significant" & music_cat != "Not_significant"]
  agree      <- sig_both[bp_cat == music_cat]

  bp_sig_set    <- xmethod[bp_cat != "Not_significant", gene]
  music_sig_set <- xmethod[music_cat != "Not_significant", gene]
  jaccard <- length(intersect(bp_sig_set, music_sig_set)) /
             length(union(bp_sig_set, music_sig_set))

  # Focus on hepatocyte-intrinsic concordance
  bp_hep    <- xmethod[bp_cat == "Hepatocyte_intrinsic", gene]
  music_hep <- xmethod[music_cat == "Hepatocyte_intrinsic", gene]
  jaccard_hep <- length(intersect(bp_hep, music_hep)) /
                 length(union(bp_hep, music_hep))

  # Build summary metrics table for dot plot
  # Also load cell-type proportion correlations if available
  concordance_path <- file.path(BP_DIR, "method_concordance.csv")
  corr_metrics <- data.table(
    metric = character(), value = numeric(), label = character()
  )

  if (file.exists(concordance_path)) {
    method_conc <- fread(concordance_path)
    if ("Hepatocytes" %in% method_conc$cell_type) {
      hep_r <- method_conc[cell_type == "Hepatocytes", pearson_r]
      hep_rho <- method_conc[cell_type == "Hepatocytes", spearman_rho]
      corr_metrics <- rbind(corr_metrics, data.table(
        metric = c("Hepatocyte\nproportion r",
                   "Hepatocyte\nproportion rho"),
        value  = c(hep_r, hep_rho),
        label  = c(sprintf("%.3f", hep_r), sprintf("%.3f", hep_rho))
      ))
    }
  }

  corr_metrics <- rbind(corr_metrics, data.table(
    metric = c("Jaccard\n(all sig)", "Jaccard\n(hep-intrinsic)",
               "Class agreement\n(both sig)"),
    value  = c(jaccard, jaccard_hep,
               if (nrow(sig_both) > 0) nrow(agree) / nrow(sig_both) else 0),
    label  = c(sprintf("%.3f", jaccard), sprintf("%.3f", jaccard_hep),
               if (nrow(sig_both) > 0) sprintf("%.1f%%", nrow(agree) / nrow(sig_both) * 100) else "N/A")
  ))

  # Order metrics for display
  corr_metrics[, metric := factor(metric, levels = rev(metric))]

  p_c <- ggplot(corr_metrics, aes(x = value, y = metric)) +
    geom_segment(aes(x = 0, xend = value, y = metric, yend = metric),
                 linewidth = 0.4, color = "gray60") +
    geom_point(size = 3, color = masld_colors$hep_intrinsic, shape = 16) +
    geom_text(aes(label = label), hjust = -0.3, size = 2.2, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    labs(x = "Concordance score",
         y = NULL,
         title = "BayesPrism vs MuSiC") +
    theme_masld() +
    theme(plot.title = element_text(size = 7, face = "bold"))
}

# ==========================================================================
# Panel d: Top 20 hepatocyte-intrinsic genes (lollipop plot)
# ==========================================================================
if (has_cols) {
  hep_genes <- deconv[category == "Hepatocyte_intrinsic" & !grepl("^ENS(G|MUSG)", symbol)]

  if (nrow(hep_genes) > 0) {
    # Top 20 by absolute adjusted logFC
    top20 <- hep_genes[order(-abs(logFC_adj))][1:min(20, .N)]
    top20[, direction := fifelse(logFC_adj > 0, "Up", "Down")]
    # Order for display (largest at top)
    top20[, symbol := factor(symbol, levels = rev(top20$symbol))]

    dir_colors <- c(Up = masld_colors$up, Down = masld_colors$down)

    p_d <- ggplot(top20, aes(x = logFC_adj, y = symbol, color = direction)) +
      geom_segment(aes(x = 0, xend = logFC_adj, y = symbol, yend = symbol),
                   linewidth = 0.5) +
      geom_point(size = 2, shape = 16) +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "gray30") +
      scale_color_manual(values = dir_colors, name = "Direction") +
      labs(x = expression("Composition-adjusted log"[2]*"FC"),
           y = NULL,
           title = "Top 20 hepatocyte-intrinsic DEGs") +
      theme_masld() +
      theme(
        axis.text.y  = element_text(face = "italic", size = 6),
        plot.title   = element_text(size = 7, face = "bold"),
        legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = "white", color = NA, linewidth = 0)
      )
  }
}

# ==========================================================================
# Assemble: 2 rows x 2 columns (A|B / C|D)
# ==========================================================================
row1 <- p_a + p_b + plot_layout(widths = c(2, 1.5))
row2 <- p_c + p_d + plot_layout(widths = c(1.5, 2))

fig3 <- (row1 / row2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig3, OUT, width = fig_full_width, height = 7)
message("Fig 3 (Deconvolution Attribution) saved to ", OUT)
