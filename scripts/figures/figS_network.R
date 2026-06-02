#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Bayesian Multiplex Gene Network
# Panels A (global community map) and D (example gene neighborhoods)
#
# SBATCH: --partition=cpu, --mem=32G, --cpus-per-task=4, --time=48:00:00,
#         --job-name=figS_network
#
# Panel A: Global Community Map
#   - Genes as points at ForceAtlas2 coordinates, colored by macro-community
#   - Community labels at centroids, convex hulls around each community
#
# Panel D: Example Gene Neighborhoods (4 subplots: THRB, HNF4A, PNPLA3, CFLAR)
#   - igraph from pre-computed gene graph JSONs
#   - Edges colored by modality layer, nodes sized by composite degree
# ==========================================================================

t0 <- proc.time()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
  library(ggforce)
  library(igraph)
  library(ggraph)
  library(jsonlite)
})

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
# Panel A: Global Community Map
# ===========================================================================
p_a <- placeholder("Panel A: global_layout.csv not found")

layout_path <- file.path(NETWORK_DIR, "global_layout.csv")
comm_path   <- file.path(NETWORK_DIR, "communities", "community_labels.csv")

if (file.exists(layout_path) && file.exists(comm_path)) {
  layout_dt <- fread(layout_path)
  comm_dt   <- fread(comm_path)

  # Merge layout with community labels
  if (all(c("gene", "x", "y") %in% names(layout_dt)) &&
      all(c("gene", "community") %in% names(comm_dt))) {

    map_dt <- merge(layout_dt, comm_dt, by = "gene", all.x = TRUE)
    map_dt[is.na(community), community := "Unassigned"]

    # Convert community to factor for coloring
    comm_counts <- map_dt[, .N, by = community][order(-N)]
    map_dt[, community := factor(community, levels = comm_counts$community)]

    # Compute centroids for labels
    centroids <- map_dt[, .(cx = median(x), cy = median(y)), by = community]
    if ("label" %in% names(comm_dt)) {
      label_dt <- unique(comm_dt[, .(community, label)])
      centroids <- merge(centroids, label_dt, by = "community", all.x = TRUE)
      centroids[is.na(label), label := as.character(community)]
    } else {
      centroids[, label := as.character(community)]
    }

    # Build a discrete palette (up to 20 communities; recycle if more)
    n_comm <- nlevels(map_dt$community)
    comm_pal <- scales::hue_pal()(min(n_comm, 20))
    if (n_comm > 20) comm_pal <- rep(comm_pal, length.out = n_comm)
    names(comm_pal) <- levels(map_dt$community)

    p_a <- ggplot(map_dt, aes(x = x, y = y, color = community)) +
      rasterize_layer(
        geom_point(size = 0.15, alpha = 0.4, shape = 16)
      ) +
      geom_mark_hull(aes(fill = community, group = community),
                     concavity = 3, alpha = 0.08, linewidth = 0.2,
                     show.legend = FALSE) +
      geom_label_repel(data = centroids,
                       aes(x = cx, y = cy, label = label),
                       size = 1.8, fontface = "bold",
                       label.padding = 0.12, box.padding = 0.3,
                       segment.size = 0.15, max.overlaps = 30,
                       inherit.aes = FALSE, fill = "white", alpha = 0.85) +
      scale_color_manual(values = comm_pal, guide = "none") +
      scale_fill_manual(values = comm_pal, guide = "none") +
      labs(title = "Global community map (ForceAtlas2 layout)",
           x = "FA2 dimension 1", y = "FA2 dimension 2") +
      theme_masld() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.line = element_blank())
  }
}

# ===========================================================================
# Panel D: Example Gene Neighborhoods (4 subplots)
# ===========================================================================
example_genes <- c("THRB", "HNF4A", "PNPLA3", "CFLAR")
neighborhood_panels <- list()

for (gene_sym in example_genes) {
  json_path <- file.path(NETWORK_DIR, "gene_graphs", paste0(gene_sym, ".json"))

  if (!file.exists(json_path)) {
    neighborhood_panels[[gene_sym]] <- placeholder(paste0(gene_sym, "\n(data not found)"))
    next
  }

  gene_data <- fromJSON(json_path, simplifyDataFrame = TRUE)

  # Build igraph from neighbors_composite
  if (is.null(gene_data$neighbors_composite) || length(gene_data$neighbors_composite) == 0) {
    neighborhood_panels[[gene_sym]] <- placeholder(paste0(gene_sym, "\n(no neighbors)"))
    next
  }

  neighbors <- as.data.table(gene_data$neighbors_composite)

  # Expect columns: source, target, layer, posterior (or weight)
  edge_cols <- intersect(c("source", "target"), names(neighbors))
  if (length(edge_cols) < 2) {
    # Try from/to naming
    if (all(c("from", "to") %in% names(neighbors))) {
      setnames(neighbors, c("from", "to"), c("source", "target"))
    } else {
      neighborhood_panels[[gene_sym]] <- placeholder(paste0(gene_sym, "\n(unexpected format)"))
      next
    }
  }

  g <- graph_from_data_frame(neighbors[, .(source, target)], directed = FALSE)

  # Edge attribute: layer type
  if ("layer" %in% names(neighbors)) {
    E(g)$layer <- neighbors$layer
  } else {
    E(g)$layer <- "unknown"
  }

  # Node attributes: composite degree and center highlight
  all_nodes <- V(g)$name
  V(g)$is_center <- all_nodes == gene_sym
  if ("composite_degree" %in% names(gene_data)) {
    deg_map <- gene_data$composite_degree
    V(g)$comp_degree <- ifelse(all_nodes %in% names(deg_map),
                               as.numeric(deg_map[all_nodes]), degree(g))
  } else {
    V(g)$comp_degree <- degree(g)
  }

  # Map layer colors
  edge_layer_vals <- E(g)$layer
  edge_cols_mapped <- ifelse(edge_layer_vals %in% names(layer_colors),
                             layer_colors[edge_layer_vals], "#BDBDBD")

  tryCatch({
    p_sub <- ggraph(g, layout = "fr") +
      geom_edge_link(aes(color = layer), alpha = 0.5, width = 0.3,
                     show.legend = FALSE) +
      geom_node_point(aes(size = comp_degree, fill = is_center),
                      shape = 21, stroke = 0.3, color = "gray30") +
      geom_node_text(aes(label = ifelse(is_center | comp_degree >= quantile(comp_degree, 0.8),
                                        name, "")),
                     size = 1.6, repel = TRUE, max.overlaps = 15) +
      scale_edge_color_manual(values = layer_colors, na.value = "#BDBDBD") +
      scale_fill_manual(values = c("TRUE" = "#FFD600", "FALSE" = "white"),
                        guide = "none") +
      scale_size_continuous(range = c(1.5, 5), guide = "none") +
      labs(title = gene_sym) +
      theme_masld() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.line = element_blank(), axis.title = element_blank(),
            plot.title = element_text(face = "bold.italic", size = 7, hjust = 0.5))

    neighborhood_panels[[gene_sym]] <- p_sub
  }, error = function(e) {
    neighborhood_panels[[gene_sym]] <<- placeholder(paste0(gene_sym, "\n(plot error)"))
  })
}

# Assemble 2x2 grid for Panel D
p_d <- wrap_plots(neighborhood_panels, ncol = 2) +
  plot_annotation(title = "Example gene neighborhoods",
                  theme = theme(plot.title = element_text(size = 8, face = "bold")))

# ===========================================================================
# Save outputs
# ===========================================================================
out_a <- file.path(FIGS_NET_DIR, "figS_network_community_map.pdf")
out_d <- file.path(FIGS_NET_DIR, "figS_network_gene_neighborhoods.pdf")

save_fig(p_a, out_a, width = fig_full_width, height = 5.5)
message("Panel A saved to ", out_a)

save_fig(p_d, out_d, width = fig_full_width, height = 5.5)
message("Panel D saved to ", out_d)

elapsed <- (proc.time() - t0)["elapsed"]
message(sprintf("figS_network.R completed in %.1f seconds", elapsed))
