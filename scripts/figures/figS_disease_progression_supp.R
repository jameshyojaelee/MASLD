#!/usr/bin/env Rscript
# figS_disease_progression_supp.R
# Supplementary Figure: Disease Progression Extended Panels
# 6 panels (a-f): transition heatmap, ordinal trend profiles,
# NAFL vs NASH volcano, early vs late fibrosis, NAS component decomposition,
# pathway enrichment by fibrosis stage
#
# Legend: docs/manuscript/05_figure_legends.md

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS02_DIR, "figS02_disease_progression_detail.pdf")

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.5

# --- Initialize all panels with placeholders ---
p_a <- placeholder("Panel a: Stage-transition heatmap")
p_b <- placeholder("Panel b: Ordinal fibrosis gene trends")
p_c <- placeholder("Panel c: NAFL vs NASH volcano")
p_d <- placeholder("Panel d: Early vs late fibrosis")
p_e <- placeholder("Panel e: NAS component decomposition")
p_f <- placeholder("Panel f: Pathway enrichment by fibrosis stage")

# ============================================================
# Panel (a): Fibrosis Stage-Transition DEG Count Heatmap
# ============================================================
fib_consec <- load_fibrosis_consecutive()
nas_consec <- load_nas_consecutive()

if (!is.null(fib_consec) && nrow(fib_consec) > 0) {
  padj_col <- intersect(c("padj", "adj.P.Val"), names(fib_consec))[1]
  fib_sig <- fib_consec[get(padj_col) < 0.05]
  trans_counts <- fib_sig[, .N, by = contrast]
  trans_counts[, direction := "All"]

  if (nrow(trans_counts) > 0) {
    p_a <- ggplot(trans_counts, aes(x = contrast, y = N)) +
      geom_col(fill = masld_colors$mash, width = 0.7) +
      geom_text(aes(label = N), vjust = -0.3, size = GEOM_TEXT_6PT) +
      labs(x = "Fibrosis transition", y = "Number of DEGs (padj < 0.05)") +
      theme_masld() +
      theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
    cat("Panel a: Stage-transition heatmap — done\n")
  }
}

# ============================================================
# Panel (b): Ordinal Fibrosis Gene Trend Profiles (Heatmap)
# ============================================================
fib_el <- load_fibrosis_early_late()

if (!is.null(fib_el) && nrow(fib_el) > 0) {
  stage_cols <- c("F0_mean", "F1_mean", "F2_mean", "F3_mean", "F4_mean")
  stage_cols <- intersect(stage_cols, names(fib_el))

  if (length(stage_cols) >= 4) {
    fib_el_complete <- fib_el[complete.cases(fib_el[, ..stage_cols])]
    fib_el_complete[, traj_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = stage_cols]
    top_genes <- fib_el_complete[order(-traj_var)][1:min(80, nrow(fib_el_complete))]

    traj_mat <- as.matrix(top_genes[, ..stage_cols])
    rownames(traj_mat) <- top_genes$symbol

    traj_scaled <- t(scale(t(traj_mat)))
    traj_scaled[is.na(traj_scaled)] <- 0

    km <- kmeans(traj_scaled, centers = 4, nstart = 25, iter.max = 100)

    cluster_centers <- km$centers
    cluster_labels <- character(4)
    for (k in 1:4) {
      center <- cluster_centers[k, ]
      if (center[length(center)] > center[1] && all(diff(center) >= -0.2)) {
        cluster_labels[k] <- "Progressive"
      } else if (center[length(center) - 1] > 0.5 && center[1] < 0.2 && center[2] < 0.2) {
        cluster_labels[k] <- "Late-onset"
      } else if (center[2] > center[length(center)]) {
        cluster_labels[k] <- "Early-transient"
      } else {
        cluster_labels[k] <- "Regressive"
      }
    }
    if (any(duplicated(cluster_labels))) {
      cluster_labels <- paste0(cluster_labels, " (", 1:4, ")")
    }

    top_genes[, cluster := cluster_labels[km$cluster]]

    cluster_order <- c("Progressive", "Late-onset", "Early-transient", "Regressive")
    cluster_order <- intersect(cluster_order, unique(top_genes$cluster))
    if (length(cluster_order) == 0) cluster_order <- unique(top_genes$cluster)

    top_genes[, cluster_f := factor(cluster, levels = cluster_order)]
    top_genes <- top_genes[order(cluster_f, -traj_var)]

    hm_mat <- as.matrix(top_genes[, ..stage_cols])
    rownames(hm_mat) <- top_genes$symbol
    colnames(hm_mat) <- gsub("_mean", "", colnames(hm_mat))

    cluster_colors <- c(
      "Progressive" = "#C2185B", "Late-onset" = "#7B1FA2",
      "Early-transient" = "#1565C0", "Regressive" = "#00695C"
    )
    avail_colors <- cluster_colors[names(cluster_colors) %in% unique(top_genes$cluster)]
    if (length(avail_colors) < length(unique(top_genes$cluster))) {
      extra <- setdiff(unique(top_genes$cluster), names(avail_colors))
      extra_cols <- setNames(rep("#757575", length(extra)), extra)
      avail_colors <- c(avail_colors, extra_cols)
    }

    row_anno <- rowAnnotation(
      Trajectory = top_genes$cluster,
      col = list(Trajectory = avail_colors),
      annotation_name_gp = gpar(fontsize = 6),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 6, fontface = "plain"),
        labels_gp = gpar(fontsize = 6)
      )
    )

    col_fun <- colorRamp2(
      c(min(hm_mat, na.rm = TRUE), 0, max(hm_mat, na.rm = TRUE)),
      c("#1565C0", "white", "#C2185B")
    )

    hm <- Heatmap(
      hm_mat,
      name = "Mean\nExpr",
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      left_annotation = row_anno,
      row_split = factor(top_genes$cluster, levels = cluster_order),
      row_gap = unit(1, "mm"),
      row_title_gp = gpar(fontsize = 6, fontface = "plain"),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 6, fontface = "plain"),
        labels_gp = gpar(fontsize = 6),
        legend_height = unit(2, "cm")
      ),
      width = unit(3, "cm"),
      use_raster = TRUE,
      raster_quality = 3
    )

    p_b <- wrap_elements(full = grid.grabExpr(draw(hm, merge_legend = TRUE)))
    message("[caption] Fibrosis stage gene trajectories")
    cat("Panel b: Ordinal fibrosis gene trends — done\n")
  }
}

# ============================================================
# Panel (c): NAFL vs NASH Volcano with NAS-Component Coloring
# ============================================================
mash_masl <- load_mash_vs_masl_results()
nas_comp_f <- load_nas_components()

if (!is.null(mash_masl) && nrow(mash_masl) > 0) {
  mash_masl[, neg_log10p := pmin(-log10(padj), 50)]
  mash_masl[, component_driver := "No data"]

  if (!is.null(nas_comp_f) && nrow(nas_comp_f) > 0) {
    padj_col_f <- intersect(c("padj", "adj.P.Val"), names(nas_comp_f))[1]
    transition_contrasts <- c("Steatosis_vs_Normal", "Inflammation_vs_Steatosis",
                               "Fibrosis_vs_Inflammation")
    comp_t <- nas_comp_f[contrast %in% transition_contrasts,
      .(gene, contrast, abs_t = abs(t))]
    dominant <- comp_t[, .SD[which.max(abs_t)], by = gene]
    driver_map <- c(
      Steatosis_vs_Normal = "Steatosis",
      Inflammation_vs_Steatosis = "Inflammation",
      Fibrosis_vs_Inflammation = "Fibrosis"
    )
    dominant[, driver := driver_map[contrast]]
    mash_masl[dominant, component_driver := i.driver, on = "gene"]
  }

  mash_masl[, sig := padj < PADJ_THRESH & abs(logFC) > LFC_THRESH]
  mash_masl[, color_group := fifelse(!sig, "NS", component_driver)]
  color_order <- c("Steatosis", "Inflammation", "Fibrosis", "No data", "NS")
  mash_masl[, color_group := factor(color_group, levels = color_order)]

  volcano_colors <- c(
    Steatosis    = "#F57F17",
    Inflammation = "#C2185B",
    Fibrosis     = "#AD1457",
    "No data"    = "#9c9c9c",
    NS           = "#cccccc"
  )

  top_label <- mash_masl[sig == TRUE][order(padj)][1:min(20, sum(mash_masl$sig, na.rm = TRUE))]

  p_c <- ggplot(mash_masl[order(-as.integer(color_group))],
                aes(x = logFC, y = neg_log10p, color = color_group)) +
    rasterize_layer(geom_point(size = 0.3, alpha = 0.5, shape = 16)) +
    geom_hline(yintercept = -log10(PADJ_THRESH), linetype = "dashed",
               linewidth = 0.3, color = "gray50") +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
               linewidth = 0.3, color = "gray50") +
    scale_color_manual(values = volcano_colors, name = "Component\nDriver",
                       drop = FALSE) +
    labs(x = expression(log[2]~"fold change (NASH vs NAFL)"),
         y = expression(-log[10](p[adj]))) +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.size = unit(0.25, "cm"),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 6))

  if (nrow(top_label) > 0 && "symbol" %in% names(top_label)) {
    p_c <- p_c +
      geom_text_repel(
        data = top_label,
        aes(label = symbol),
        size = GEOM_TEXT_6PT,
        max.overlaps = 15,
        segment.size = 0.2,
        color = "black",
        fontface = "italic",
        min.segment.length = 0.2
      )
  }
  cat("Panel c: NAFL vs NASH volcano — done\n")
}

# ============================================================
# Panel (d): Early vs Late Fibrosis Programs
# ============================================================
if (!is.null(fib_el) && nrow(fib_el) > 0 && "stage_class" %in% names(fib_el)) {
  # Compute late-early difference from stage means
  early_cols <- intersect(c("F0_mean", "F1_mean", "F2_mean"), names(fib_el))
  late_cols  <- intersect(c("F3_mean", "F4_mean"), names(fib_el))

  if (length(early_cols) > 0 && length(late_cols) > 0) {
    fib_el_complete <- fib_el[complete.cases(fib_el[, c(..early_cols, ..late_cols)])]
    fib_el_complete[, early_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = early_cols]
    fib_el_complete[, late_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = late_cols]
    fib_el_complete[, diff_late_early := late_mean - early_mean]

    # Top 15 Early + top 15 Late by absolute difference
    top_early <- fib_el_complete[stage_class == "Early"][order(diff_late_early)][1:min(15, sum(fib_el_complete$stage_class == "Early"))]
    top_late  <- fib_el_complete[stage_class == "Late"][order(-diff_late_early)][1:min(15, sum(fib_el_complete$stage_class == "Late"))]
    top_el <- rbindlist(list(top_early, top_late))
    top_el[, direction := fifelse(stage_class == "Late", "Late-enriched", "Early-enriched")]
    top_el <- top_el[order(diff_late_early)]
    top_el[, symbol := factor(symbol, levels = symbol)]

    p_d <- ggplot(top_el, aes(x = diff_late_early, y = symbol, color = direction)) +
      geom_point(size = 1.5, shape = 16) +
      geom_segment(aes(x = 0, xend = diff_late_early, y = symbol, yend = symbol),
                   linewidth = 0.3) +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "gray50") +
      scale_color_manual(values = c("Early-enriched" = masld_colors$control,
                                     "Late-enriched" = masld_colors$fibrosis),
                         name = NULL) +
      labs(x = "Mean expression diff (Late - Early)", y = NULL) +
      theme_masld() +
      theme(axis.text.y = element_text(size = 6, face = "italic"),
            legend.position = "bottom",
            legend.key.size = unit(0.25, "cm"))

    cat("Panel d: Early vs late fibrosis — done\n")
  }
}

# ============================================================
# Panel (e): NAS Component Decomposition
# ============================================================
nas_comp <- load_nas_components()

if (!is.null(nas_comp) && nrow(nas_comp) > 0) {
  padj_col <- intersect(c("padj", "adj.P.Val"), names(nas_comp))[1]
  nas_sig <- nas_comp[get(padj_col) < PADJ_THRESH & abs(logFC) > LFC_THRESH]

  comp_up <- nas_sig[logFC > 0, .(n_degs = .N, direction = "Up"), by = contrast]
  comp_down <- nas_sig[logFC < 0, .(n_degs = -.N, direction = "Down"), by = contrast]
  comp_bar <- rbindlist(list(comp_up, comp_down))

  contrast_order <- c("Steatosis_vs_Normal", "Inflammation_vs_Steatosis",
                       "Fibrosis_vs_Inflammation", "Ordinal_progression")
  contrast_labels <- c("Steatosis\nvs Normal", "Inflammation\nvs Steatosis",
                        "Fibrosis\nvs Inflammation", "Ordinal\nProgression")

  comp_bar <- comp_bar[contrast %in% contrast_order]
  comp_bar[, contrast_f := factor(contrast, levels = contrast_order,
                                   labels = contrast_labels)]

  component_colors <- c(
    "Steatosis\nvs Normal"          = "#F57F17",
    "Inflammation\nvs Steatosis"    = "#C2185B",
    "Fibrosis\nvs Inflammation"     = "#AD1457",
    "Ordinal\nProgression"          = "#7B1FA2"
  )

  p_e <- ggplot(comp_bar, aes(x = contrast_f, y = n_degs, fill = contrast_f)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = component_colors) +
    labs(x = "Disease Transition", y = "Number of DEGs") +
    theme_masld() +
    theme(axis.text.x = element_text(size = 6, lineheight = 0.9))

  cat("Panel e: NAS component decomposition — done\n")
}

# ============================================================
# Panel (f): Pathway Switching Across Fibrosis Stages
# ============================================================
fib_gsea <- load_fibrosis_stage_gsea()

if (!is.null(fib_gsea) && nrow(fib_gsea) > 0) {
  sig_pathways <- fib_gsea[padj < 0.1, unique(pathway)]

  if (length(sig_pathways) > 0) {
    pathway_rank <- fib_gsea[pathway %in% sig_pathways,
      .(min_padj = min(padj, na.rm = TRUE)), by = pathway][order(min_padj)]
    top_pathways <- pathway_rank$pathway[1:min(20, nrow(pathway_rank))]

    gsea_plot <- fib_gsea[pathway %in% top_pathways]
    gsea_plot[, pathway_clean := gsub("HALLMARK_", "", pathway)]
    gsea_plot[, pathway_clean := gsub("_", " ", pathway_clean)]
    gsea_plot[, pathway_clean := tools::toTitleCase(tolower(pathway_clean))]

    gsea_plot[, stage_label := gsub("_vs_F0", "", contrast)]
    gsea_plot[, stage_label := factor(stage_label, levels = c("F1", "F2", "F3", "F4"))]
    gsea_plot[, neg_log10p := pmin(-log10(padj), 10)]

    pathway_peak <- gsea_plot[, .(peak_stage = stage_label[which.max(abs(NES))]),
                              by = pathway_clean]
    pathway_peak[, peak_num := as.integer(peak_stage)]
    pathway_order <- pathway_peak[order(peak_num, decreasing = TRUE), pathway_clean]
    gsea_plot[, pathway_clean := factor(pathway_clean, levels = pathway_order)]

    p_f <- ggplot(gsea_plot, aes(x = stage_label, y = pathway_clean)) +
      geom_point(aes(size = neg_log10p, color = NES), shape = 16) +
      scale_size_continuous(
        name = expression(-log[10](padj)),
        range = c(0.5, 4),
        breaks = c(1, 2, 5, 10)
      ) +
      scale_color_gradient2(
        low = "#1565C0", mid = "white", high = "#C2185B",
        midpoint = 0, name = "NES",
        limits = c(-max(abs(gsea_plot$NES)), max(abs(gsea_plot$NES)))
      ) +
      labs(x = "Fibrosis Stage (vs F0)", y = NULL) +
      theme_masld() +
      theme(
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6)
      )

    cat("Panel f: Pathway switching dot plot — done\n")
  }
}

# ============================================================
# ASSEMBLE FIGURE
# ============================================================
cat("\nAssembling supplementary figure 8...\n")

# Layout: 3 rows of 2
fig <- (p_a | p_b) /
       (p_c | p_d) /
       (p_e | p_f) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 6, face = "plain"))

save_fig_tall(fig, OUT, width = fig_full_width, height = 9.5, dpi = 300)
cat("Supplementary Figure (disease progression) saved to:", OUT, "\n")
