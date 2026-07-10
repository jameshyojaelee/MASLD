#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Bayesian Multiplex Network — Statistics
# Panels B (informativeness heatmap), C (multiplicity), F (diagnostics)
#
# SBATCH: --partition=cpu, --mem=16G, --cpus-per-task=4, --time=48:00:00,
#         --job-name=figS_network_stats
#
# Panel B: Modality Informativeness Heatmap
#   - Mean posterior per macro-community per layer (ComplexHeatmap)
#
# Panel C: Edge Multiplicity Distribution
#   - Stacked bar: x = K (1-10), y = count of gene pairs, color = max layer
#
# Panel F: Bayesian Diagnostics
#   - 10-panel grid: null vs observed density per modality layer
# ==========================================================================

t0 <- proc.time()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
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
# Panel B: Modality Informativeness Heatmap
# ===========================================================================
out_b <- file.path(FIGS_NET_DIR, "figS_network_informativeness.pdf")

# Look for posterior edge data + community assignments
posterior_path <- file.path(NETWORK_DIR, "posterior_edges.csv")
comm_path      <- file.path(NETWORK_DIR, "communities", "community_labels.csv")

if (file.exists(posterior_path) && file.exists(comm_path)) {
  edges_dt <- fread(posterior_path)
  comm_dt  <- fread(comm_path)

  # Require gene_a, gene_b, layer, posterior columns
  needed_edge <- c("gene_a", "gene_b", "layer", "posterior")
  # Try alternate column names
  if (!"gene_a" %in% names(edges_dt) && "source" %in% names(edges_dt))
    setnames(edges_dt, "source", "gene_a")
  if (!"gene_b" %in% names(edges_dt) && "target" %in% names(edges_dt))
    setnames(edges_dt, "target", "gene_b")
  if (!"posterior" %in% names(edges_dt) && "weight" %in% names(edges_dt))
    setnames(edges_dt, "weight", "posterior")

  if (all(needed_edge %in% names(edges_dt)) &&
      all(c("gene", "community") %in% names(comm_dt))) {

    # Assign community to each gene in edge
    comm_map <- comm_dt[, .(gene, community)]
    edges_dt <- merge(edges_dt, comm_map, by.x = "gene_a", by.y = "gene", all.x = TRUE)
    setnames(edges_dt, "community", "comm_a")
    edges_dt <- merge(edges_dt, comm_map, by.x = "gene_b", by.y = "gene", all.x = TRUE)
    setnames(edges_dt, "community", "comm_b")

    # Keep intra-community edges only (community-level informativeness)
    intra <- edges_dt[!is.na(comm_a) & !is.na(comm_b) & comm_a == comm_b]

    if (nrow(intra) > 0) {
      # Mean posterior per community x layer
      info_mat <- dcast(intra, comm_a ~ layer, value.var = "posterior",
                        fun.aggregate = mean, fill = 0)
      rn <- info_mat$comm_a
      info_mat[, comm_a := NULL]
      mat <- as.matrix(info_mat)
      rownames(mat) <- rn

      # Use community labels if available
      if ("label" %in% names(comm_dt)) {
        label_map <- unique(comm_dt[, .(community, label)])
        row_labels <- label_map$label[match(rownames(mat), label_map$community)]
        row_labels[is.na(row_labels)] <- rownames(mat)[is.na(row_labels)]
        rownames(mat) <- row_labels
      }

      # Color function
      col_fun <- colorRamp2(c(0, 0.5, 1), c("#E3F2FD", "#42A5F5", "#0D47A1"))

      # Column annotation with layer colors
      col_layers <- colnames(mat)
      col_anno_colors <- ifelse(col_layers %in% names(layer_colors),
                                layer_colors[col_layers], "#BDBDBD")
      ha <- HeatmapAnnotation(
        Layer = col_layers,
        col = list(Layer = setNames(col_anno_colors, col_layers)),
        show_legend = FALSE,
        annotation_name_gp = gpar(fontsize = 6)
      )

      ht <- Heatmap(mat,
                    name = "Mean\nposterior",
                    col = col_fun,
                    top_annotation = ha,
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    column_names_rot = 45,
                    row_title = "Community",
                    column_title_gp = gpar(fontsize = 6, fontface = "plain"),
                    heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                labels_gp = gpar(fontsize = 6)))

      message("[caption] Modality layer informativeness")
      pdf(out_b, width = fig_full_width, height = 5)
      draw(ht, padding = unit(c(2, 2, 2, 10), "mm"))
      dev.off()
      message("Panel B saved to ", out_b)
    } else {
      message("Panel B: no intra-community edges found")
      pdf(out_b, width = fig_full_width, height = 5)
      grid.newpage()
      grid.text("Panel B: no intra-community edges found", gp = gpar(fontsize = 6))
      dev.off()
    }
  } else {
    message("Panel B: required columns not found in posterior_edges.csv")
    pdf(out_b, width = fig_full_width, height = 5)
    grid.newpage()
    grid.text("Panel B: required columns not found", gp = gpar(fontsize = 6))
    dev.off()
  }
} else {
  message("Panel B: posterior_edges.csv or community_labels.csv not found")
  pdf(out_b, width = fig_full_width, height = 5)
  grid.newpage()
  grid.text("Panel B: data files not found", gp = gpar(fontsize = 6))
  dev.off()
}

# ===========================================================================
# Panel C: Edge Multiplicity Distribution
# ===========================================================================
p_c <- placeholder("Panel C: composite edges not found")

composite_path <- file.path(NETWORK_DIR, "composite_edges.csv")
if (file.exists(composite_path)) {
  comp_dt <- fread(composite_path)

  # Expect k_multiplicity column and a max-layer column
  if (!"k_multiplicity" %in% names(comp_dt)) {
    # Try to derive from per-layer posterior columns
    layer_cols <- intersect(names(layer_colors), names(comp_dt))
    if (length(layer_cols) > 0) {
      comp_dt[, k_multiplicity := rowSums(.SD > 0.5, na.rm = TRUE), .SDcols = layer_cols]
      comp_dt[, max_layer := layer_cols[max.col(.SD, ties.method = "first")], .SDcols = layer_cols]
    }
  }

  if ("k_multiplicity" %in% names(comp_dt)) {
    # Determine max contributing layer if not present
    if (!"max_layer" %in% names(comp_dt)) {
      layer_cols <- intersect(names(layer_colors), names(comp_dt))
      if (length(layer_cols) > 0) {
        comp_dt[, max_layer := layer_cols[max.col(.SD, ties.method = "first")], .SDcols = layer_cols]
      } else {
        comp_dt[, max_layer := "unknown"]
      }
    }

    # Summarize by K x max_layer
    mult_summary <- comp_dt[, .N, by = .(k_multiplicity, max_layer)]
    mult_summary[, k_multiplicity := factor(k_multiplicity)]

    # Use layer colors
    used_layers <- unique(mult_summary$max_layer)
    bar_colors <- ifelse(used_layers %in% names(layer_colors),
                         layer_colors[used_layers], "#BDBDBD")
    names(bar_colors) <- used_layers

    p_c <- ggplot(mult_summary, aes(x = k_multiplicity, y = N, fill = max_layer)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = bar_colors, name = "Max-contributing\nlayer") +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
      labs(x = "Edge multiplicity (K)",
           y = "Number of gene pairs") +
      theme_masld() +
      theme(legend.position = "right")
  }
}

message("[caption] Edge multiplicity distribution")
out_c <- file.path(FIGS_NET_DIR, "figS_network_multiplicity.pdf")
save_fig(p_c, out_c, width = fig_full_width, height = 3.5)
message("Panel C saved to ", out_c)

# ===========================================================================
# Panel F: Bayesian Diagnostics — null vs observed density per layer
# ===========================================================================
p_f <- placeholder("Panel F: bayesian_diagnostics.csv not found")

diag_path <- file.path(NETWORK_DIR, "bayesian_diagnostics.csv")
null_dir  <- file.path(NETWORK_DIR, "null_distributions")

if (file.exists(diag_path) && dir.exists(null_dir)) {
  diag_dt <- fread(diag_path)

  # Load null distribution density files
  null_files <- list.files(null_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(null_files) > 0) {
    null_list <- lapply(null_files, function(f) {
      dt <- fread(f)
      # Extract layer name from filename
      layer_name <- gsub("\\.csv$", "", basename(f))
      layer_name <- gsub("^null_", "", layer_name)
      dt[, layer := layer_name]
      dt[, distribution := "null"]
      dt
    })
    null_dt <- rbindlist(null_list, fill = TRUE)

    # Expect columns: x (or value), density, layer
    if (!"x" %in% names(null_dt) && "value" %in% names(null_dt))
      setnames(null_dt, "value", "x")

    # Build observed density from diagnostics
    # diag_dt should have layer + posterior columns
    if ("layer" %in% names(diag_dt) && "posterior" %in% names(diag_dt) &&
        all(c("x", "density") %in% names(null_dt))) {

      # Compute observed kernel density per layer
      obs_list <- lapply(split(diag_dt, by = "layer"), function(sub) {
        if (nrow(sub) < 5) return(NULL)
        d <- density(sub$posterior, from = 0, to = 1, n = 256)
        data.table(x = d$x, density = d$y, layer = sub$layer[1],
                   distribution = "observed")
      })
      obs_dt <- rbindlist(obs_list, fill = TRUE)

      # Combine null and observed
      plot_dt <- rbind(null_dt[, .(x, density, layer, distribution)],
                       obs_dt, fill = TRUE)

      # Extract pi_0 from diagnostics if available
      pi0_dt <- NULL
      if ("pi_0" %in% names(diag_dt)) {
        pi0_dt <- unique(diag_dt[, .(layer, pi_0)])
      }

      # Keep only layers present in both
      shared_layers <- intersect(unique(null_dt$layer), unique(obs_dt$layer))
      plot_dt <- plot_dt[layer %in% shared_layers]

      # Color map
      dist_colors <- c(null = "#BDBDBD", observed = "#C2185B")

      p_f <- ggplot(plot_dt, aes(x = x, y = density, color = distribution, fill = distribution)) +
        geom_area(data = plot_dt[distribution == "null"],
                  alpha = 0.3, linewidth = 0.3) +
        geom_line(data = plot_dt[distribution == "observed"],
                  linewidth = 0.5) +
        scale_color_manual(values = dist_colors, name = NULL) +
        scale_fill_manual(values = dist_colors, name = NULL, guide = "none") +
        facet_wrap(~ layer, ncol = 5, scales = "free_y") +
        labs(x = "Edge weight", y = "Density") +
        theme_masld() +
        theme(legend.position = "bottom",
              strip.text = element_text(size = 6, face = "plain"))

      # Add pi_0 annotations if available
      if (!is.null(pi0_dt) && nrow(pi0_dt) > 0) {
        pi0_dt <- pi0_dt[layer %in% shared_layers]
        if (nrow(pi0_dt) > 0) {
          # Get y range for annotation placement
          p_f <- p_f +
            geom_text(data = pi0_dt,
                      aes(x = 0.75, y = Inf, label = sprintf("pi[0]==%.2f", pi_0)),
                      parse = TRUE, inherit.aes = FALSE,
                      vjust = 1.5, hjust = 0.5, size = GEOM_TEXT_6PT, color = "black")
        }
      }
    }
  }
}

message("[caption] Bayesian diagnostics: null vs observed edge weight distributions")
out_f <- file.path(FIGS_NET_DIR, "figS_network_diagnostics.pdf")
save_fig(p_f, out_f, width = fig_full_width, height = 5)
message("Panel F saved to ", out_f)

elapsed <- (proc.time() - t0)["elapsed"]
message(sprintf("figS_network_stats.R completed in %.1f seconds", elapsed))
