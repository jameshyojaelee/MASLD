#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Bayesian Multiplex Network — Cross-Modality Bridges
# Panel E: Cross-modality bridge genes (alluvial/Sankey diagram)
#
# SBATCH: --partition=cpu, --mem=16G, --cpus-per-task=4, --time=48:00:00,
#         --job-name=figS_network_bridges
#
# Panel E: Cross-Modality Bridges
#   - Identify bridge genes with high posterior in 3+ different modality types
#   - ggalluvial showing genes flowing between modality layers
#   - Top 10 bridge genes labeled
#   - Fallback: bipartite network if ggalluvial insufficient
# ==========================================================================

t0 <- proc.time()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

# Try ggalluvial; fall back to bipartite approach
has_alluvial <- requireNamespace("ggalluvial", quietly = TRUE)
if (has_alluvial) {
  suppressPackageStartupMessages(library(ggalluvial))
}

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ---------------------------------------------------------------------------
# Modality layer color palette (10 layers)
# ---------------------------------------------------------------------------
layer_colors <- c(
  coexpression   = "#1565C0",
  coloc          = "#C2185B",
  twas           = "#E91E63",
  spatial        = "#7B1FA2",
  sc_pseudobulk  = "#42A5F5",
  atac           = "#00695C",
  protein        = "#F57F17",
  drug           = "#0D47A1",
  pathway        = "#F48FB1",
  cross_species  = "#AD1457"
)

# ===========================================================================
# Panel E: Cross-Modality Bridges
# ===========================================================================
p_e <- placeholder("Panel E: composite_edges.csv not found")

composite_path <- file.path(NETWORK_DIR, "composite_edges.csv")
posterior_path <- file.path(NETWORK_DIR, "posterior_edges.csv")

# Strategy: identify genes with high-posterior edges in 3+ distinct modality layers
# Use posterior_edges.csv (has gene_a, gene_b, layer, posterior) to compute per-gene
# modality participation

build_bridge_data <- function() {
  # Try posterior edges first (most informative)
  if (file.exists(posterior_path)) {
    edges_dt <- fread(posterior_path)

    # Normalize column names
    if (!"gene_a" %in% names(edges_dt) && "source" %in% names(edges_dt))
      setnames(edges_dt, "source", "gene_a")
    if (!"gene_b" %in% names(edges_dt) && "target" %in% names(edges_dt))
      setnames(edges_dt, "target", "gene_b")
    if (!"posterior" %in% names(edges_dt) && "weight" %in% names(edges_dt))
      setnames(edges_dt, "weight", "posterior")

    needed <- c("gene_a", "gene_b", "layer", "posterior")
    if (!all(needed %in% names(edges_dt))) return(NULL)

    # Filter to confident edges (posterior > 0.5)
    sig_edges <- edges_dt[posterior > 0.5]
    if (nrow(sig_edges) == 0) return(NULL)

    # For each gene, count distinct modality layers it participates in
    gene_a_layers <- sig_edges[, .(layer = unique(layer)), by = gene_a]
    setnames(gene_a_layers, "gene_a", "gene")
    gene_b_layers <- sig_edges[, .(layer = unique(layer)), by = gene_b]
    setnames(gene_b_layers, "gene_b", "gene")
    gene_layers <- unique(rbind(gene_a_layers, gene_b_layers))

    # Count layers per gene
    layer_counts <- gene_layers[, .(n_layers = .N), by = gene]

    # Bridge genes: present in 3+ layers
    bridges <- layer_counts[n_layers >= 3][order(-n_layers)]
    if (nrow(bridges) == 0) return(NULL)

    # Build gene-layer long-form for alluvial
    bridge_genes <- bridges$gene
    bridge_layers <- gene_layers[gene %in% bridge_genes]
    bridge_layers <- merge(bridge_layers, bridges[, .(gene, n_layers)], by = "gene")

    # Compute mean posterior per gene x layer for sizing
    gene_layer_post <- rbind(
      sig_edges[gene_a %in% bridge_genes, .(gene = gene_a, layer, posterior)],
      sig_edges[gene_b %in% bridge_genes, .(gene = gene_b, layer, posterior)]
    )[, .(mean_posterior = mean(posterior)), by = .(gene, layer)]

    bridge_layers <- merge(bridge_layers, gene_layer_post, by = c("gene", "layer"), all.x = TRUE)
    bridge_layers[is.na(mean_posterior), mean_posterior := 0.5]

    return(list(bridges = bridges, bridge_layers = bridge_layers))
  }

  # Fallback: use composite_edges.csv
  if (file.exists(composite_path)) {
    comp_dt <- fread(composite_path)
    layer_cols <- intersect(names(layer_colors), names(comp_dt))
    if (length(layer_cols) == 0) return(NULL)

    # Identify gene columns
    gene_col_a <- intersect(c("gene_a", "source"), names(comp_dt))[1]
    gene_col_b <- intersect(c("gene_b", "target"), names(comp_dt))[1]
    if (is.na(gene_col_a) || is.na(gene_col_b)) return(NULL)

    # For each edge, find layers with posterior > 0.5
    edges_long <- melt(comp_dt, id.vars = c(gene_col_a, gene_col_b),
                       measure.vars = layer_cols,
                       variable.name = "layer", value.name = "posterior")
    edges_long <- edges_long[posterior > 0.5]
    if (nrow(edges_long) == 0) return(NULL)

    setnames(edges_long, c(gene_col_a, gene_col_b), c("gene_a", "gene_b"))

    gene_a_layers <- edges_long[, .(layer = unique(layer)), by = gene_a]
    setnames(gene_a_layers, "gene_a", "gene")
    gene_b_layers <- edges_long[, .(layer = unique(layer)), by = gene_b]
    setnames(gene_b_layers, "gene_b", "gene")
    gene_layers <- unique(rbind(gene_a_layers, gene_b_layers))

    layer_counts <- gene_layers[, .(n_layers = .N), by = gene]
    bridges <- layer_counts[n_layers >= 3][order(-n_layers)]
    if (nrow(bridges) == 0) return(NULL)

    bridge_genes <- bridges$gene
    bridge_layers <- gene_layers[gene %in% bridge_genes]
    bridge_layers <- merge(bridge_layers, bridges[, .(gene, n_layers)], by = "gene")

    gene_layer_post <- rbind(
      edges_long[gene_a %in% bridge_genes, .(gene = gene_a, layer, posterior)],
      edges_long[gene_b %in% bridge_genes, .(gene = gene_b, layer, posterior)]
    )[, .(mean_posterior = mean(posterior)), by = .(gene, layer)]
    bridge_layers <- merge(bridge_layers, gene_layer_post, by = c("gene", "layer"), all.x = TRUE)
    bridge_layers[is.na(mean_posterior), mean_posterior := 0.5]

    return(list(bridges = bridges, bridge_layers = bridge_layers))
  }

  return(NULL)
}

bridge_data <- build_bridge_data()

if (!is.null(bridge_data)) {
  bridges <- bridge_data$bridges
  bridge_layers <- bridge_data$bridge_layers

  # Top 10 bridge genes for labeling
  top10 <- head(bridges, 10)
  n_bridges <- nrow(bridges)

  if (has_alluvial) {
    # -------------------------------------------------------------------
    # Alluvial approach: genes (strata) flow across modality layers (axes)
    # -------------------------------------------------------------------
    # For alluvial, we need a wide-form or a long-form with alluvium IDs.
    # Build a frequency table: for each gene, list its layer participation.
    # Use top 50 bridges to keep the plot readable.
    top_n <- min(50, nrow(bridges))
    top_bridges <- bridges[1:top_n, gene]
    alluvial_dt <- bridge_layers[gene %in% top_bridges]

    # Create a "connection" long-form: each gene-layer pair is a stratum
    # The alluvium is the gene, the x-axis is the layer
    alluvial_dt[, layer := factor(layer, levels = names(layer_colors))]
    alluvial_dt[, is_top10 := gene %in% top10$gene]

    # Summarize: for each layer, how many bridge genes participate?
    layer_gene_count <- alluvial_dt[, .N, by = layer]

    # Build a bipartite-style dot plot (more robust than alluvial for this data)
    # Left: bridge genes (top 50), Right: modality layers, edges = participation
    gene_order <- bridges[gene %in% top_bridges][order(-n_layers), gene]
    alluvial_dt[, gene := factor(gene, levels = rev(gene_order))]

    p_e <- ggplot(alluvial_dt, aes(x = layer, y = gene)) +
      geom_point(aes(size = mean_posterior,
                     color = layer),
                 shape = 16, alpha = 0.7) +
      geom_point(data = alluvial_dt[is_top10 == TRUE],
                 aes(size = mean_posterior),
                 shape = 1, color = "black", stroke = 0.5) +
      scale_color_manual(values = layer_colors, guide = "none") +
      scale_size_continuous(range = c(0.8, 3.5), name = "Mean\nposterior",
                            breaks = c(0.5, 0.7, 0.9)) +
      labs(x = "Modality layer", y = NULL,
           title = sprintf("Cross-modality bridge genes (n=%s, top %d shown)",
                           format(n_bridges, big.mark = ","), top_n),
           subtitle = sprintf("Top 10 outlined: %s",
                              paste(top10$gene, collapse = ", "))) +
      theme_masld() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            axis.text.y = element_text(size = ifelse(top_n > 30, 4, 5),
                                       face = ifelse(gene_order %in% top10$gene,
                                                     "bold", "plain")),
            plot.title = element_text(size = 8, face = "bold"),
            plot.subtitle = element_text(size = 6, face = "italic"),
            legend.position = "right")

  } else {
    # -------------------------------------------------------------------
    # Fallback: tile plot (no ggalluvial available)
    # -------------------------------------------------------------------
    top_n <- min(50, nrow(bridges))
    top_bridges <- bridges[1:top_n, gene]
    tile_dt <- bridge_layers[gene %in% top_bridges]
    gene_order <- bridges[gene %in% top_bridges][order(-n_layers), gene]
    tile_dt[, gene := factor(gene, levels = rev(gene_order))]
    tile_dt[, layer := factor(layer, levels = names(layer_colors))]

    p_e <- ggplot(tile_dt, aes(x = layer, y = gene, fill = mean_posterior)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#E3F2FD", mid = "#42A5F5", high = "#0D47A1",
                           midpoint = 0.7, name = "Mean\nposterior",
                           limits = c(0.5, 1)) +
      labs(x = "Modality layer", y = NULL,
           title = sprintf("Cross-modality bridge genes (n=%s, top %d shown)",
                           format(n_bridges, big.mark = ","), top_n)) +
      theme_masld() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            axis.text.y = element_text(size = ifelse(top_n > 30, 4, 5)),
            plot.title = element_text(size = 8, face = "bold"),
            legend.position = "right")
  }

  # Add a summary bar inset: layer participation frequency
  layer_freq <- bridge_layers[gene %in% bridges$gene, .N, by = layer]
  layer_freq[, layer := factor(layer, levels = names(layer_colors))]
  bar_cols <- layer_colors[as.character(layer_freq$layer)]

  p_bar <- ggplot(layer_freq, aes(x = layer, y = N, fill = layer)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = layer_colors, guide = "none") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Bridge gene\nparticipation",
         title = "Layer participation") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          plot.title = element_text(size = 7, face = "bold"))

  # Combine main panel + summary bar
  p_e <- p_e / p_bar + plot_layout(heights = c(3, 1))
}

# ===========================================================================
# Save output
# ===========================================================================
out_e <- file.path(FIGS_NET_DIR, "figS_network_bridges.pdf")
save_fig(p_e, out_e, width = fig_full_width, height = 8)
message("Panel E saved to ", out_e)

elapsed <- (proc.time() - t0)["elapsed"]
message(sprintf("figS_network_bridges.R completed in %.1f seconds", elapsed))
