#!/usr/bin/env Rscript
#' 307: Pseudotime Publication Figures.
#'
#' Generates 7 supplementary PDF figure panels for the scRNA-seq pseudotime
#' analysis pipeline (Scripts 300-306). Gracefully skips any panel whose
#' input file is missing.
#'
#' Inputs (from results_gpu_v2/pseudotime/):
#'   - consensus_pseudotime_all.csv  (306)
#'   - cytotrace2_scores_{ct}.csv    (302)
#'   - robust_trajectory_genes.csv   (306)
#'   - gene_dynamics_{ct}.csv        (301)
#'   - stage_switch_genes_{ct}.csv   (301)
#'   - paga_connectivity_all.csv     (303)
#'   - paga_connectivity_healthy.csv (303)
#'   - paga_connectivity_nash.csv    (303)
#'   - paga_cluster_annotations.csv  (303)
#'   - method_concordance.csv        (306)
#'   - palantir_pseudotime_{ct}.csv  (301)
#'   - monocle3_pseudotime_{ct}.csv  (304)
#'   - bulk_sc_concordance.csv       (305)
#'
#' Outputs (to figures/supplementary/figS02_progression/):
#'   - figS_pseudotime_umap.pdf
#'   - figS_pseudotime_potency.pdf
#'   - figS_pseudotime_gene_heatmap.pdf
#'   - figS_pseudotime_switch_genes.pdf
#'   - figS_paga_connectivity.pdf
#'   - figS_pseudotime_concordance.pdf
#'   - figS_bulk_sc_validation.pdf
#'
#' Usage:
#'   Rscript 307_pseudotime_figures.R
#'   # or via SLURM: sbatch run_pseudotime_pipeline.sh

# =========================================================================
# Libraries
# =========================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridis)
  library(igraph)
})

# ComplexHeatmap is Bioconductor — load with fallback
has_complex_heatmap <- suppressPackageStartupMessages(
  requireNamespace("ComplexHeatmap", quietly = TRUE)
)
if (has_complex_heatmap) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
  })
} else {
  message("ComplexHeatmap not installed — heatmap panels will use ggplot2 fallback")
}

# =========================================================================
# Paths
# =========================================================================
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

PT_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime")
FIG_DIR <- file.path(BASE, "figures/supplementary/figS02_progression")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                "Endothelial_cells", "Cholangiocytes")

# Display-friendly names for facet labels
CT_LABELS <- c(
  Hepatocytes       = "Hepatocytes",
  Macrophages       = "Macrophages",
  Fibroblasts       = "Fibroblasts",
  Endothelial_cells = "Endothelial cells",
  Cholangiocytes    = "Cholangiocytes"
)

cat("==========================================================\n")
cat("307: Pseudotime Publication Figures\n")
cat("PT_DIR: ", PT_DIR, "\n")
cat("FIG_DIR:", FIG_DIR, "\n")
cat("==========================================================\n")


# =========================================================================
# Helper: safe file loader
# =========================================================================
safe_fread <- function(path, ...) {
  if (!file.exists(path)) {
    warning("File not found: ", path, call. = FALSE)
    return(NULL)
  }
  tryCatch(fread(path, ...), error = function(e) {
    warning("Failed to read: ", path, " — ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

# =========================================================================
# Condition ordering used throughout
# =========================================================================
CONDITION_LEVELS <- c("Healthy", "NAFLD", "MASLD", "NASH", "Cirrhotic")

condition_colors <- c(
  Healthy   = masld_colors$control,
  NAFLD     = masld_colors$nafl,
  MASLD     = "#F06292",
  NASH      = masld_colors$nash,
  Cirrhotic = masld_colors$fibrosis
)

# Cell-type colors for PAGA and multi-panel figures
ct_fill_colors <- c(
  Hepatocytes       = "#0D47A1",
  Macrophages       = "#C2185B",
  Fibroblasts       = "#F57F17",
  Endothelial_cells = "#2E7D32",
  Cholangiocytes    = "#1565C0"
)


# =========================================================================
# PANEL 1: UMAP colored by consensus pseudotime
# =========================================================================
cat("\n--- Panel 1: Pseudotime UMAP ---\n")

consensus_path <- file.path(PT_DIR, "consensus_pseudotime_all.csv")
consensus_df   <- safe_fread(consensus_path)

if (!is.null(consensus_df)) {

  # Identify UMAP columns (may be UMAP_1/UMAP_2 or umap_1/umap_2)
  umap_cols <- grep("^[Uu][Mm][Aa][Pp]", names(consensus_df), value = TRUE)
  pt_col    <- grep("consensus.*pseudotime|pseudotime.*consensus",
                    names(consensus_df), value = TRUE, ignore.case = TRUE)
  ct_col    <- grep("cell_type|celltype", names(consensus_df),
                    value = TRUE, ignore.case = TRUE)

  # Fallback: if no UMAP columns in CSV, try loading from h5ad-exported coords
  if (length(umap_cols) < 2) {
    # Attempt to construct UMAP from per-celltype palantir/dpt outputs
    umap_pieces <- list()
    for (ct in CELL_TYPES) {
      # Try Palantir output (has cell index, can merge with subset UMAP)
      pal_path <- file.path(PT_DIR, paste0("palantir_pseudotime_", ct, ".csv"))
      if (file.exists(pal_path)) {
        pal_df <- fread(pal_path)
        # First column is cell barcode (index)
        if ("V1" %in% names(pal_df)) setnames(pal_df, "V1", "cell")
        pal_df[, cell_type := ct]
        umap_pieces[[ct]] <- pal_df[, .(cell, cell_type)]
      }
    }
    if (length(umap_pieces) > 0) {
      message("  UMAP coordinates not found in consensus CSV — ",
              "UMAP panel will use pseudotime only (no spatial layout)")
    }
  }

  # Find pseudotime column robustly
  if (length(pt_col) == 0) {
    pt_col <- grep("pseudotime", names(consensus_df), value = TRUE,
                   ignore.case = TRUE)
  }
  if (length(ct_col) == 0) {
    ct_col <- grep("cell.type", names(consensus_df), value = TRUE,
                   ignore.case = TRUE)
  }

  can_plot_umap <- length(umap_cols) >= 2 && length(pt_col) >= 1 && length(ct_col) >= 1

  if (can_plot_umap) {
    plot_df <- consensus_df[, c(umap_cols[1:2], pt_col[1], ct_col[1]), with = FALSE]
    setnames(plot_df, c("UMAP1", "UMAP2", "pseudotime", "cell_type"))

    # Clean cell type labels
    plot_df[, cell_type_label := CT_LABELS[cell_type]]
    plot_df <- plot_df[!is.na(cell_type_label)]

    # Order by pseudotime so high values are plotted on top
    plot_df <- plot_df[order(pseudotime)]

    p_umap <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
      geom_point(size = 0.05, alpha = 0.6) +
      scale_color_viridis(option = "viridis", name = "Consensus\npseudotime") +
      facet_wrap(~cell_type_label, nrow = 1, scales = "free") +
      labs(x = "UMAP 1", y = "UMAP 2",
           title = "Consensus pseudotime across cell types") +
      theme_masld(base_size = 7) +
      theme(legend.position = "right",
            legend.key.width = unit(0.2, "cm"),
            legend.key.height = unit(0.5, "cm"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(size = 7, face = "bold"))

    outfile <- file.path(FIG_DIR, "figS_pseudotime_umap.pdf")
    save_fig(p_umap, outfile, width = fig_full_width, height = 3)
    cat("  Saved:", outfile, "\n")

  } else {
    message("  Skipping UMAP panel: insufficient columns (need UMAP + pseudotime + cell_type)")
  }

} else {
  message("  Skipping Panel 1: consensus_pseudotime_all.csv not found")
}


# =========================================================================
# PANEL 2: CytoTRACE potency violins by condition and cell type
# =========================================================================
cat("\n--- Panel 2: CytoTRACE potency violins ---\n")

potency_pieces <- list()
for (ct in CELL_TYPES) {
  path <- file.path(PT_DIR, paste0("cytotrace2_scores_", ct, ".csv"))
  if (!file.exists(path)) path <- file.path(PT_DIR, paste0("cytotrace_v1_scores_", ct, ".csv"))
  df <- safe_fread(path)
  if (!is.null(df)) {
    # First column may be cell index
    if ("V1" %in% names(df)) setnames(df, "V1", "cell")
    df[, cell_type := ct]
    potency_pieces[[ct]] <- df
  }
}

if (length(potency_pieces) > 0) {
  potency_all <- rbindlist(potency_pieces, fill = TRUE)

  # Filter to known conditions and set factor levels
  potency_all <- potency_all[condition %in% CONDITION_LEVELS]
  potency_all[, condition := factor(condition, levels = CONDITION_LEVELS)]
  potency_all[, cell_type_label := CT_LABELS[cell_type]]

  # Kruskal-Wallis per cell type
  kw_results <- potency_all[, {
    grps <- split(potency, condition)
    grps <- grps[sapply(grps, length) >= 5]
    if (length(grps) >= 2) {
      kw <- kruskal.test(grps)
      .(kw_p = kw$p.value)
    } else {
      .(kw_p = NA_real_)
    }
  }, by = cell_type]
  kw_results[, p_label := fifelse(
    is.na(kw_p), "",
    fifelse(kw_p < 0.001, "***",
    fifelse(kw_p < 0.01, "**",
    fifelse(kw_p < 0.05, "*", "ns")))
  )]
  kw_results[, cell_type_label := CT_LABELS[cell_type]]

  # Compute y position for annotation
  y_max <- potency_all[, .(ymax = max(potency, na.rm = TRUE) * 1.05),
                        by = cell_type_label]
  kw_results <- merge(kw_results, y_max, by = "cell_type_label", all.x = TRUE)

  p_potency <- ggplot(potency_all, aes(x = condition, y = potency, fill = condition)) +
    geom_violin(linewidth = 0.2, scale = "width", alpha = 0.7, show.legend = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.2,
                 fill = "white", alpha = 0.6, show.legend = FALSE) +
    geom_text(data = kw_results,
              aes(x = 2.5, y = ymax, label = paste0("KW: ", p_label)),
              inherit.aes = FALSE, size = 2, hjust = 0.5, vjust = -0.2) +
    scale_fill_manual(values = condition_colors) +
    facet_wrap(~cell_type_label, nrow = 1, scales = "free_y") +
    labs(x = "Disease condition", y = "CytoTRACE potency",
         title = "Cell potency by disease condition and cell type") +
    theme_masld(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          strip.text = element_text(size = 7, face = "bold"))

  outfile <- file.path(FIG_DIR, "figS_pseudotime_potency.pdf")
  save_fig(p_potency, outfile, width = fig_full_width, height = 3.5)
  cat("  Saved:", outfile, "\n")

} else {
  message("  Skipping Panel 2: no cytotrace2_scores files found")
}


# =========================================================================
# PANEL 3: Gene dynamics heatmap — top 50 trajectory genes per cell type
# =========================================================================
cat("\n--- Panel 3: Gene dynamics heatmap ---\n")

robust_genes_df <- safe_fread(file.path(PT_DIR, "robust_trajectory_genes.csv"))

if (!is.null(robust_genes_df)) {

  # Identify gene and cell_type columns
  gene_col <- grep("^gene$", names(robust_genes_df), value = TRUE, ignore.case = TRUE)
  ct_col_r <- grep("cell_type|celltype", names(robust_genes_df), value = TRUE,
                   ignore.case = TRUE)

  if (length(gene_col) == 0) gene_col <- names(robust_genes_df)[1]
  if (length(ct_col_r) == 0) ct_col_r <- grep("cell", names(robust_genes_df),
                                                value = TRUE, ignore.case = TRUE)[1]

  # Select top 50 per cell type (by absolute correlation or rank)
  rank_col <- grep("rank|rho|corr", names(robust_genes_df), value = TRUE,
                   ignore.case = TRUE)

  heatmap_list <- list()

  for (ct in CELL_TYPES) {
    dyn_path <- file.path(PT_DIR, paste0("gene_dynamics_", ct, ".csv"))
    dyn_df   <- safe_fread(dyn_path)
    if (is.null(dyn_df)) next

    # gene_dynamics has gene as first col (index), then bin_0 .. bin_N
    gene_names <- dyn_df[[1]]
    bin_cols <- grep("^bin_", names(dyn_df), value = TRUE)
    if (length(bin_cols) == 0) next

    mat <- as.matrix(dyn_df[, ..bin_cols])
    rownames(mat) <- gene_names

    # Select top 50 from robust_trajectory_genes for this cell type
    if (!is.null(ct_col_r) && !is.na(ct_col_r)) {
      ct_genes <- robust_genes_df[get(ct_col_r) == ct, get(gene_col)]
    } else {
      ct_genes <- robust_genes_df[[gene_col]]
    }
    ct_genes <- ct_genes[ct_genes %in% gene_names]
    ct_genes <- head(ct_genes, 50)

    if (length(ct_genes) < 5) next

    mat_sub <- mat[ct_genes, , drop = FALSE]

    # Z-score rows for visualization
    mat_z <- t(scale(t(mat_sub)))
    mat_z[is.na(mat_z)] <- 0

    heatmap_list[[ct]] <- mat_z
  }

  if (length(heatmap_list) > 0 && has_complex_heatmap) {
    outfile <- file.path(FIG_DIR, "figS_pseudotime_gene_heatmap.pdf")
    pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
    pdf_device(outfile, width = fig_full_width, height = 2.5 * length(heatmap_list))

    ht_list <- NULL
    for (ct in names(heatmap_list)) {
      ht <- Heatmap(
        heatmap_list[[ct]],
        name = paste0(CT_LABELS[ct], "\nz-score"),
        col = colorRamp2(c(-2, 0, 2), c(masld_colors$down, "white", masld_colors$up)),
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 4, fontface = "italic"),
        column_title = CT_LABELS[ct],
        column_title_gp = gpar(fontsize = 7, fontface = "bold"),
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 6),
          labels_gp = gpar(fontsize = 5),
          legend_width = unit(2, "cm"),
          direction = "horizontal"
        ),
        row_title = paste0("Top ", nrow(heatmap_list[[ct]]), " genes"),
        row_title_gp = gpar(fontsize = 6),
        height = unit(2, "cm") * (nrow(heatmap_list[[ct]]) / 50),
        use_raster = TRUE, raster_quality = 3
      )
      if (is.null(ht_list)) {
        ht_list <- ht
      } else {
        ht_list <- ht_list %v% ht
      }
    }

    # Add pseudotime axis label
    column_ha <- HeatmapAnnotation(
      pseudotime = anno_simple(
        seq(0, 1, length.out = ncol(heatmap_list[[1]])),
        col = colorRamp2(c(0, 1), c("#FDE725", "#440154"))
      ),
      annotation_name_gp = gpar(fontsize = 6),
      show_annotation_name = TRUE,
      annotation_label = "Pseudotime",
      height = unit(3, "mm")
    )

    draw(ht_list, heatmap_legend_side = "bottom",
         column_title = "Top trajectory genes ordered by pseudotime",
         column_title_gp = gpar(fontsize = 8, fontface = "bold"))

    dev.off()
    cat("  Saved:", outfile, "\n")

  } else if (length(heatmap_list) > 0) {
    # ggplot2 tile fallback if ComplexHeatmap unavailable
    tile_pieces <- list()
    for (ct in names(heatmap_list)) {
      mat_z <- heatmap_list[[ct]]
      tile_df <- as.data.frame(mat_z) %>%
        tibble::rownames_to_column("gene") %>%
        pivot_longer(-gene, names_to = "bin", values_to = "zscore") %>%
        mutate(cell_type = ct,
               bin_num = as.integer(gsub("bin_", "", bin)))
      tile_pieces <- c(tile_pieces, list(tile_df))
    }
    tile_all <- bind_rows(tile_pieces)
    tile_all$cell_type_label <- CT_LABELS[tile_all$cell_type]

    p_heat <- ggplot(tile_all, aes(x = bin_num, y = gene, fill = zscore)) +
      geom_tile() +
      scale_fill_gradient2(low = masld_colors$down, mid = "white",
                           high = masld_colors$up, midpoint = 0,
                           name = "Z-score", limits = c(-2, 2),
                           oob = scales::squish) +
      facet_wrap(~cell_type_label, scales = "free_y", ncol = 1) +
      labs(x = "Pseudotime bin", y = NULL,
           title = "Top trajectory genes ordered by pseudotime") +
      theme_masld(base_size = 7) +
      theme(axis.text.y = element_text(size = 3, face = "italic"),
            strip.text = element_text(size = 7, face = "bold"),
            legend.key.height = unit(0.3, "cm"))

    outfile <- file.path(FIG_DIR, "figS_pseudotime_gene_heatmap.pdf")
    save_fig_tall(p_heat, outfile, width = fig_full_width,
                  height = 2.5 * length(heatmap_list))
    cat("  Saved:", outfile, "\n")
  }

} else {
  message("  Skipping Panel 3: robust_trajectory_genes.csv not found")
}


# =========================================================================
# PANEL 4: Switch genes heatmap per cell type
# =========================================================================
cat("\n--- Panel 4: Switch genes heatmap ---\n")

switch_pieces <- list()
for (ct in CELL_TYPES) {
  path <- file.path(PT_DIR, paste0("stage_switch_genes_", ct, ".csv"))
  df <- safe_fread(path)
  if (!is.null(df) && nrow(df) > 0) {
    df[, cell_type := ct]
    switch_pieces[[ct]] <- df
  }
}

if (length(switch_pieces) > 0) {
  switch_all <- rbindlist(switch_pieces, fill = TRUE)

  # Top 20 switch genes per cell type (by adjusted p-value)
  switch_top <- switch_all[, head(.SD[order(padj)], 20), by = cell_type]

  # Identify required columns
  has_lfc <- "lfc" %in% names(switch_top)
  has_bin <- "bin" %in% names(switch_top)
  has_cond <- "bin_condition" %in% names(switch_top)

  if (has_lfc && has_bin) {
    switch_top[, cell_type_label := CT_LABELS[cell_type]]
    switch_top[, neg_log10_padj := -log10(pmax(padj, 1e-300))]

    # Create a tile heatmap: genes x bins, colored by LFC
    p_switch <- ggplot(switch_top,
                       aes(x = factor(bin), y = reorder(gene, -padj),
                           fill = lfc)) +
      geom_tile(color = "white", linewidth = 0.1) +
      scale_fill_gradient2(low = masld_colors$down, mid = "white",
                           high = masld_colors$up, midpoint = 0,
                           name = "Log FC") +
      facet_wrap(~cell_type_label, scales = "free", ncol = 3) +
      labs(x = "Pseudotime bin", y = NULL,
           title = "Disease-stage switch genes") +
      theme_masld(base_size = 7) +
      theme(axis.text.y = element_text(size = 4, face = "italic"),
            axis.text.x = element_text(size = 5),
            strip.text = element_text(size = 7, face = "bold"),
            legend.key.height = unit(0.3, "cm"))

    # Add condition annotation bar if available
    if (has_cond) {
      p_switch <- p_switch +
        labs(subtitle = "Bins labeled by dominant disease condition")
    }

    n_ct <- length(unique(switch_top$cell_type))
    fig_height <- max(4, 2.5 * ceiling(n_ct / 3))

    outfile <- file.path(FIG_DIR, "figS_pseudotime_switch_genes.pdf")
    save_fig(p_switch, outfile, width = fig_full_width, height = fig_height)
    cat("  Saved:", outfile, "\n")

  } else {
    message("  Skipping Panel 4: missing expected columns (lfc, bin)")
  }

} else {
  message("  Skipping Panel 4: no stage_switch_genes files found")
}


# =========================================================================
# PANEL 5: PAGA connectivity graph
# =========================================================================
cat("\n--- Panel 5: PAGA connectivity ---\n")

paga_all_path     <- file.path(PT_DIR, "paga_connectivity_all.csv")
paga_healthy_path <- file.path(PT_DIR, "paga_connectivity_healthy.csv")
paga_nash_path    <- file.path(PT_DIR, "paga_connectivity_nash.csv")
paga_annot_path   <- file.path(PT_DIR, "paga_cluster_annotations.csv")

paga_all_df   <- safe_fread(paga_all_path)
paga_annot_df <- safe_fread(paga_annot_path)

if (!is.null(paga_all_df) && nrow(paga_all_df) > 0 && !is.null(paga_annot_df)) {

  # --- Helper to build an igraph from PAGA edge list ---
  build_paga_graph <- function(edge_df, annot_df, min_weight = 0.1) {
    edge_df <- edge_df[weight >= min_weight]
    if (nrow(edge_df) == 0) return(NULL)

    g <- graph_from_data_frame(
      edge_df[, .(source, target, weight)],
      directed = FALSE,
      vertices = annot_df[, .(name = leiden_fine)]
    )

    # Assign cell-type colors
    ct_map <- setNames(annot_df$dominant_celltype, annot_df$leiden_fine)
    V(g)$celltype <- ct_map[V(g)$name]
    V(g)$color <- ct_fill_colors[V(g)$celltype]
    V(g)$color[is.na(V(g)$color)] <- "gray70"

    # Node size proportional to cell count
    n_map <- setNames(annot_df$n_cells, annot_df$leiden_fine)
    V(g)$size <- sqrt(as.numeric(n_map[V(g)$name])) / 5
    V(g)$size[is.na(V(g)$size)] <- 2

    # Edge width proportional to weight
    E(g)$width <- E(g)$weight * 3

    return(g)
  }

  # Ensure leiden_fine is character
  paga_annot_df[, leiden_fine := as.character(leiden_fine)]
  paga_all_df[, source := as.character(source)]
  paga_all_df[, target := as.character(target)]

  g_all <- build_paga_graph(paga_all_df, paga_annot_df)

  # Build Healthy and NASH graphs
  paga_healthy_df <- safe_fread(paga_healthy_path)
  paga_nash_df    <- safe_fread(paga_nash_path)

  if (!is.null(paga_healthy_df) && nrow(paga_healthy_df) > 0) {
    paga_healthy_df[, source := as.character(source)]
    paga_healthy_df[, target := as.character(target)]
  }
  if (!is.null(paga_nash_df) && nrow(paga_nash_df) > 0) {
    paga_nash_df[, source := as.character(source)]
    paga_nash_df[, target := as.character(target)]
  }

  g_healthy <- if (!is.null(paga_healthy_df) && nrow(paga_healthy_df) > 0) {
    build_paga_graph(paga_healthy_df, paga_annot_df)
  } else NULL

  g_nash <- if (!is.null(paga_nash_df) && nrow(paga_nash_df) > 0) {
    build_paga_graph(paga_nash_df, paga_annot_df)
  } else NULL

  # Determine layout: how many panels
  n_panels <- 1 + !is.null(g_healthy) + !is.null(g_nash)

  outfile <- file.path(FIG_DIR, "figS_paga_connectivity.pdf")
  pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
  pdf_device(outfile, width = 3.5 * n_panels, height = 4)

  par(mfrow = c(1, n_panels), mar = c(1, 1, 2, 1), family = "Helvetica")

  # Compute stable layout from the full graph
  set.seed(42)
  if (!is.null(g_all)) {
    layout_all <- layout_with_fr(g_all)

    plot(g_all, layout = layout_all,
         vertex.label = NA,
         vertex.frame.color = "gray30",
         vertex.frame.width = 0.3,
         edge.color = adjustcolor("gray40", alpha.f = 0.5),
         main = "All conditions")

    # Legend
    ct_present <- unique(na.omit(V(g_all)$celltype))
    legend("bottomleft", legend = CT_LABELS[ct_present],
           col = ct_fill_colors[ct_present], pch = 16, cex = 0.55,
           bty = "n", pt.cex = 1)
  }

  if (!is.null(g_healthy)) {
    # Use same layout coordinates (subset of nodes)
    nodes_h <- V(g_healthy)$name
    idx_h <- match(nodes_h, V(g_all)$name)
    layout_h <- layout_all[idx_h, , drop = FALSE]
    # Handle missing nodes with fresh coordinates
    missing <- is.na(idx_h)
    if (any(missing)) {
      layout_h[missing, ] <- matrix(runif(sum(missing) * 2, -5, 5),
                                     ncol = 2)
    }

    plot(g_healthy, layout = layout_h,
         vertex.label = NA,
         vertex.frame.color = "gray30",
         vertex.frame.width = 0.3,
         edge.color = adjustcolor(masld_colors$down, alpha.f = 0.5),
         main = "Healthy")
  }

  if (!is.null(g_nash)) {
    nodes_n <- V(g_nash)$name
    idx_n <- match(nodes_n, V(g_all)$name)
    layout_n <- layout_all[idx_n, , drop = FALSE]
    missing <- is.na(idx_n)
    if (any(missing)) {
      layout_n[missing, ] <- matrix(runif(sum(missing) * 2, -5, 5),
                                     ncol = 2)
    }

    plot(g_nash, layout = layout_n,
         vertex.label = NA,
         vertex.frame.color = "gray30",
         vertex.frame.width = 0.3,
         edge.color = adjustcolor(masld_colors$up, alpha.f = 0.5),
         main = "NASH")
  }

  dev.off()
  cat("  Saved:", outfile, "\n")

} else {
  message("  Skipping Panel 5: paga_connectivity files not found")
}


# =========================================================================
# PANEL 6: Cross-method concordance
# =========================================================================
cat("\n--- Panel 6: Method concordance ---\n")

concordance_path <- file.path(PT_DIR, "method_concordance.csv")
concordance_df   <- safe_fread(concordance_path)

if (!is.null(concordance_df) && nrow(concordance_df) > 0) {

  # Identify columns
  ct_col_c   <- grep("cell_type|celltype", names(concordance_df),
                     value = TRUE, ignore.case = TRUE)
  rho_col    <- grep("spearman|rho|correlation", names(concordance_df),
                     value = TRUE, ignore.case = TRUE)
  method_col <- grep("method|comparison", names(concordance_df),
                     value = TRUE, ignore.case = TRUE)

  if (length(ct_col_c) > 0 && length(rho_col) > 0) {
    plot_df <- concordance_df[, c(ct_col_c[1], rho_col[1]), with = FALSE]
    setnames(plot_df, c("cell_type", "spearman_rho"))

    if (length(method_col) > 0) {
      plot_df[, method := concordance_df[[method_col[1]]]]
    } else {
      plot_df[, method := "Cross-method"]
    }

    plot_df[, cell_type_label := CT_LABELS[cell_type]]
    plot_df <- plot_df[!is.na(cell_type_label)]

    # Panel A: Bar chart of correlations
    pA_conc <- ggplot(plot_df, aes(x = reorder(cell_type_label, spearman_rho),
                                    y = spearman_rho, fill = method)) +
      geom_col(width = 0.6, alpha = 0.85) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60",
                 linewidth = 0.3) +
      coord_flip() +
      scale_fill_brewer(palette = "Set2", name = "Comparison") +
      labs(x = NULL, y = "Spearman rho",
           title = "Cross-method pseudotime agreement") +
      theme_masld(base_size = 7) +
      theme(legend.position = "right",
            legend.key.size = unit(0.25, "cm"))
  } else {
    pA_conc <- placeholder("Concordance: missing rho column")
  }

  # Panel B: Palantir vs Monocle 3 scatter (if both available per cell type)
  scatter_panels <- list()
  for (ct in CELL_TYPES) {
    pal_path <- file.path(PT_DIR, paste0("palantir_pseudotime_", ct, ".csv"))
    m3_path  <- file.path(PT_DIR, paste0("monocle3_pseudotime_", ct, ".csv"))

    pal_df <- safe_fread(pal_path)
    m3_df  <- safe_fread(m3_path)

    if (!is.null(pal_df) && !is.null(m3_df)) {
      # Index column
      if ("V1" %in% names(pal_df)) setnames(pal_df, "V1", "cell")
      pal_pt_col <- grep("palantir.*pseudotime", names(pal_df),
                         value = TRUE, ignore.case = TRUE)

      m3_cell_col <- grep("^cell$", names(m3_df), value = TRUE, ignore.case = TRUE)
      m3_pt_col   <- grep("monocle3.*pseudotime", names(m3_df),
                          value = TRUE, ignore.case = TRUE)

      if (length(pal_pt_col) > 0 && length(m3_pt_col) > 0) {
        # Merge
        if (length(m3_cell_col) > 0) {
          setnames(m3_df, m3_cell_col[1], "cell")
        } else if ("V1" %in% names(m3_df)) {
          setnames(m3_df, "V1", "cell")
        }

        # Ensure cell column in palantir
        if (!"cell" %in% names(pal_df)) {
          pal_df[, cell := .I]  # fallback
        }

        merged <- merge(
          pal_df[, c("cell", pal_pt_col[1]), with = FALSE],
          m3_df[, c("cell", m3_pt_col[1]), with = FALSE],
          by = "cell"
        )
        setnames(merged, c("cell", "palantir", "monocle3"))
        merged <- merged[is.finite(palantir) & is.finite(monocle3)]

        if (nrow(merged) > 50) {
          rho_val <- cor(merged$palantir, merged$monocle3, method = "spearman")

          # Subsample for plotting
          set.seed(42)
          plot_sub <- merged[sample(.N, min(.N, 5000))]

          p_scat <- ggplot(plot_sub, aes(x = palantir, y = monocle3)) +
            geom_point(size = 0.2, alpha = 0.3, color = masld_colors$up) +
            geom_smooth(method = "lm", se = FALSE, linewidth = 0.4,
                        color = "black", linetype = "dashed") +
            annotate("text", x = Inf, y = -Inf,
                     label = sprintf("rho = %.3f", rho_val),
                     hjust = 1.1, vjust = -0.5, size = 2.2) +
            labs(x = "Palantir pseudotime", y = "Monocle 3 pseudotime",
                 title = CT_LABELS[ct]) +
            theme_masld(base_size = 7)

          scatter_panels[[ct]] <- p_scat
        }
      }
    }
  }

  # Assemble concordance figure
  if (length(scatter_panels) > 0) {
    scatter_combined <- wrap_plots(scatter_panels, nrow = 1)
    combined_conc <- pA_conc / scatter_combined +
      plot_layout(heights = c(1, 1)) +
      plot_annotation(tag_levels = "a") &
      theme(plot.tag = element_text(size = 8, face = "bold"))

    outfile <- file.path(FIG_DIR, "figS_pseudotime_concordance.pdf")
    n_scat <- length(scatter_panels)
    save_fig(combined_conc, outfile,
             width = max(fig_full_width, 1.5 * n_scat),
             height = 5.5)
    cat("  Saved:", outfile, "\n")

  } else {
    # Just bar chart
    outfile <- file.path(FIG_DIR, "figS_pseudotime_concordance.pdf")
    save_fig(pA_conc, outfile, width = fig_half_width, height = 3.5)
    cat("  Saved:", outfile, "\n")
  }

} else {
  # Fallback: build concordance from raw Palantir + Monocle3 + DPT files
  cat("  method_concordance.csv not found — building from raw pseudotime files\n")

  corr_rows <- list()
  for (ct in CELL_TYPES) {
    pal_path <- file.path(PT_DIR, paste0("palantir_pseudotime_", ct, ".csv"))
    dpt_path <- file.path(PT_DIR, paste0("dpt_pseudotime_", ct, ".csv"))
    m3_path  <- file.path(PT_DIR, paste0("monocle3_pseudotime_", ct, ".csv"))

    pal_df <- safe_fread(pal_path)
    dpt_df <- safe_fread(dpt_path)
    m3_df  <- safe_fread(m3_path)

    # Standardize index column
    standardize_idx <- function(df) {
      if (is.null(df)) return(NULL)
      if ("V1" %in% names(df)) setnames(df, "V1", "cell")
      if ("cell" %in% names(df)) return(df)
      # Use first column as cell ID
      setnames(df, names(df)[1], "cell")
      df
    }

    pal_df <- standardize_idx(pal_df)
    dpt_df <- standardize_idx(dpt_df)
    m3_df  <- standardize_idx(m3_df)

    # Extract pseudotime columns
    extract_pt <- function(df, pattern) {
      if (is.null(df)) return(NULL)
      col <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
      if (length(col) == 0) return(NULL)
      data.table(cell = df$cell, pt = df[[col[1]]])
    }

    pt_pal <- extract_pt(pal_df, "palantir.*pseudotime")
    pt_dpt <- extract_pt(dpt_df, "dpt.*pseudotime")
    pt_m3  <- extract_pt(m3_df, "monocle3.*pseudotime")

    pt_list <- list(Palantir = pt_pal, DPT = pt_dpt, Monocle3 = pt_m3)
    pt_list <- pt_list[!sapply(pt_list, is.null)]

    if (length(pt_list) >= 2) {
      combos <- combn(names(pt_list), 2)
      for (k in seq_len(ncol(combos))) {
        m1 <- combos[1, k]
        m2 <- combos[2, k]
        merged <- merge(pt_list[[m1]], pt_list[[m2]], by = "cell",
                        suffixes = c("_1", "_2"))
        merged <- merged[is.finite(pt_1) & is.finite(pt_2)]
        if (nrow(merged) > 30) {
          rho <- cor(merged$pt_1, merged$pt_2, method = "spearman")
          corr_rows[[length(corr_rows) + 1]] <- data.table(
            cell_type = ct,
            method = paste0(m1, " vs ", m2),
            spearman_rho = rho,
            n_cells = nrow(merged)
          )
        }
      }
    }
  }

  if (length(corr_rows) > 0) {
    conc_built <- rbindlist(corr_rows)
    conc_built[, cell_type_label := CT_LABELS[cell_type]]
    conc_built <- conc_built[!is.na(cell_type_label)]

    p_conc_fb <- ggplot(conc_built,
                        aes(x = reorder(cell_type_label, spearman_rho),
                            y = spearman_rho, fill = method)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60",
                 linewidth = 0.3) +
      coord_flip() +
      scale_fill_brewer(palette = "Set2", name = "Comparison") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = NULL, y = "Spearman rho",
           title = "Cross-method pseudotime agreement") +
      theme_masld(base_size = 7) +
      theme(legend.position = "right",
            legend.key.size = unit(0.25, "cm"))

    outfile <- file.path(FIG_DIR, "figS_pseudotime_concordance.pdf")
    save_fig(p_conc_fb, outfile, width = fig_full_width, height = 3.5)
    cat("  Saved:", outfile, "\n")
  } else {
    message("  Skipping Panel 6: no pseudotime files found for concordance")
  }
}


# =========================================================================
# PANEL 7: Bulk-SC concordance validation
# =========================================================================
cat("\n--- Panel 7: Bulk-SC validation ---\n")

bulk_sc_path <- file.path(PT_DIR, "bulk_sc_concordance.csv")
bulk_sc_df   <- safe_fread(bulk_sc_path)

if (!is.null(bulk_sc_df) && nrow(bulk_sc_df) > 0) {

  # Identify columns
  ct_col_b  <- grep("cell_type|celltype", names(bulk_sc_df),
                    value = TRUE, ignore.case = TRUE)
  sig_col   <- grep("signature|comparison", names(bulk_sc_df),
                    value = TRUE, ignore.case = TRUE)
  rho_col_b <- grep("spearman|rho", names(bulk_sc_df),
                    value = TRUE, ignore.case = TRUE)
  jac_col   <- grep("jaccard", names(bulk_sc_df),
                    value = TRUE, ignore.case = TRUE)

  panels_7 <- list()

  # Panel A: Spearman rho bar chart (signature-pseudotime correlations)
  if (length(ct_col_b) > 0 && length(rho_col_b) > 0) {
    rho_df <- bulk_sc_df[!is.na(get(rho_col_b[1]))]

    if (length(sig_col) > 0) {
      rho_df[, signature_short := gsub("^sig_|^transition_", "",
                                        get(sig_col[1]))]
    } else {
      rho_df[, signature_short := paste0("Sig_", seq_len(.N))]
    }

    rho_df[, cell_type_label := CT_LABELS[get(ct_col_b[1])]]
    rho_df <- rho_df[!is.na(cell_type_label)]

    # Cap display to top 20 per cell type if too many
    rho_df <- rho_df[, head(.SD[order(-abs(get(rho_col_b[1])))], 20),
                     by = cell_type_label]

    pA_bsc <- ggplot(rho_df, aes(x = reorder(signature_short,
                                              get(rho_col_b[1])),
                                  y = get(rho_col_b[1]),
                                  fill = get(rho_col_b[1]) > 0)) +
      geom_col(width = 0.7, show.legend = FALSE) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      coord_flip() +
      scale_fill_manual(values = c(`TRUE` = masld_colors$up,
                                    `FALSE` = masld_colors$down)) +
      facet_wrap(~cell_type_label, scales = "free_y", ncol = 3) +
      labs(x = NULL, y = "Spearman rho (signature vs pseudotime)",
           title = "Bulk signature correlation with single-cell pseudotime") +
      theme_masld(base_size = 7) +
      theme(axis.text.y = element_text(size = 4.5),
            strip.text = element_text(size = 7, face = "bold"))

    panels_7[["rho"]] <- pA_bsc
  }

  # Panel B: Jaccard overlap bar chart
  if (length(jac_col) > 0 && length(ct_col_b) > 0) {
    jac_df <- bulk_sc_df[!is.na(get(jac_col[1]))]
    jac_df[, cell_type_label := CT_LABELS[get(ct_col_b[1])]]
    jac_df <- jac_df[!is.na(cell_type_label)]

    if (nrow(jac_df) > 0) {
      pB_bsc <- ggplot(jac_df, aes(x = cell_type_label,
                                    y = get(jac_col[1]),
                                    fill = cell_type_label)) +
        geom_col(width = 0.6, show.legend = FALSE) +
        scale_fill_manual(values = ct_fill_colors) +
        labs(x = NULL, y = "Jaccard index",
             title = "Gene overlap: SC trajectory vs bulk transport") +
        theme_masld(base_size = 7) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

      panels_7[["jaccard"]] <- pB_bsc
    }
  }

  # Assemble
  if (length(panels_7) == 2) {
    combined_bsc <- panels_7[["rho"]] / panels_7[["jaccard"]] +
      plot_layout(heights = c(2, 1)) +
      plot_annotation(tag_levels = "a") &
      theme(plot.tag = element_text(size = 8, face = "bold"))

    outfile <- file.path(FIG_DIR, "figS_bulk_sc_validation.pdf")
    save_fig(combined_bsc, outfile, width = fig_full_width, height = 6)
    cat("  Saved:", outfile, "\n")

  } else if (length(panels_7) == 1) {
    outfile <- file.path(FIG_DIR, "figS_bulk_sc_validation.pdf")
    save_fig(panels_7[[1]], outfile, width = fig_full_width, height = 4)
    cat("  Saved:", outfile, "\n")
  }

} else {
  message("  Skipping Panel 7: bulk_sc_concordance.csv not found")
}


# =========================================================================
# Done
# =========================================================================
cat("\n==========================================================\n")
cat("307: Pseudotime figures COMPLETE\n")
cat("Output directory:", FIG_DIR, "\n")
cat("==========================================================\n")
