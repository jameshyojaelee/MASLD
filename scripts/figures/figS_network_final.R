#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Final Bayesian Multiplex Gene Network (4 panels)
# Post-remediation pipeline (16,704 genes, 7 active layers, 3.6M composite pairs)
#
# SBATCH: --partition=cpu, --mem=32G, --cpus-per-task=4, --time=48:00:00,
#         --job-name=figS_net_final
#
# Panel A: Global Community Map (FA2 layout, macro communities + hallmark labels)
# Panel B: Modality Informativeness Heatmap (communities x 7 layers)
# Panel C: Gold-Standard Recovery Benchmark (fold over null with 95% CI)
# Panel D: Example Gene Neighborhoods (THRB, HNF4A, PNPLA3, CFLAR)
#
# Dependencies: ID-normalized network outputs (network_nodes.csv,
# composite_edges.csv, edges_{layer}.csv, communities/, gene_graphs/*.json,
# id_normalization_summary.csv). Benchmark: figures/misc/network_benchmark_v2/.
# ==========================================================================

t0 <- proc.time()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(ggforce)
  library(igraph)
  library(ggraph)
  library(jsonlite)
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
# Paths
# ---------------------------------------------------------------------------
NET_DIR       <- NETWORK_DIR
COMM_DIR      <- file.path(NET_DIR, "communities")
GRAPH_DIR     <- file.path(NET_DIR, "gene_graphs")
BENCH_DIR     <- file.path(BASE, "figures/misc/network_benchmark_v2/data")
OUT_DIR       <- FIGS_NET_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ID normalization must be complete before this script runs
id_norm_path <- file.path(NET_DIR, "id_normalization_summary.csv")
if (!file.exists(id_norm_path)) {
  stop("id_normalization_summary.csv not found -- ID normalization pipeline must ",
       "run before this figure script. Expected: ", id_norm_path)
}

# ---------------------------------------------------------------------------
# Active layers (7): ppi, coexpr, regulon, lr, genetic, pathway, cerna
# Spatial, cosmos, xspecies layers were dropped (0 edges post-remediation)
# ---------------------------------------------------------------------------
ACTIVE_LAYERS <- c("ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna")

layer_colors <- c(
  ppi     = "#1565C0",   # Deep blue
  coexpr  = "#42A5F5",   # Light blue
  regulon = "#E91E63",   # Pink
  lr      = "#F57F17",   # Amber
  genetic = "#C2185B",   # Magenta
  pathway = "#F48FB1",   # Soft pink
  cerna   = "#7B1FA2"    # Violet
)
layer_labels <- c(
  ppi     = "PPI",
  coexpr  = "Co-expression",
  regulon = "Regulon",
  lr      = "Ligand-receptor",
  genetic = "Genetic",
  pathway = "Pathway",
  cerna   = "ceRNA"
)

# Shared community color palette
community_pal <- c(
  "0"  = "#C2185B",  # inflammation / EMT
  "1"  = "#1565C0",  # OxPhos / metabolism
  "2"  = "#7B1FA2",  # cell cycle / MYC
  "3"  = "#00695C",  # other meso
  "4"  = "#F57F17",  # other meso
  "5"  = "#AD1457",
  "6"  = "#42A5F5",
  "Other" = "#BDBDBD"
)

# Hard-coded macro labels (derived from community_enrichment.csv, hallmark,
# top-ranked terms; captured in review spec).
MACRO_LABELS <- c(
  "0" = "Inflammation / EMT / Drug targets",
  "1" = "OxPhos / Metabolism",
  "2" = "Cell cycle / MYC targets",
  "3" = "Macro 3",
  "4" = "Macro 4",
  "5" = "Macro 5",
  "6" = "Macro 6"
)

MIN_COMM_SIZE <- 30   # Hide/group very small communities

# ===========================================================================
# Load core data
# ===========================================================================
message("Loading nodes/edges/communities ...")
nodes <- fread(file.path(NET_DIR, "network_nodes.csv"))
setnames(nodes, "human_symbol", "gene")

comm <- fread(file.path(COMM_DIR, "community_assignments.csv"))
comm[, macro_id := as.character(macro_id)]
comm_sizes <- comm[, .N, by = macro_id][order(-N)]
keep_macros <- comm_sizes[N >= MIN_COMM_SIZE, macro_id]
comm[, macro_label := ifelse(macro_id %in% keep_macros, macro_id, "Other")]

enr <- fread(file.path(COMM_DIR, "community_enrichment.csv"))
enr_hallmark <- enr[category == "hallmark" & level == "macro"]
# Top hallmark term per macro community
setorder(enr_hallmark, community_id, padj)
top_hallmark <- enr_hallmark[, .SD[1], by = community_id]
top_hallmark[, community_id := as.character(community_id)]

# Merge labels: prefer MACRO_LABELS, else top hallmark, else "Macro {id}"
label_lookup <- data.table(
  macro_id = unique(comm$macro_label)
)
label_lookup[, label := ifelse(
  macro_id %in% names(MACRO_LABELS),
  MACRO_LABELS[macro_id],
  paste0("Macro ", macro_id)
)]
# Fill in hallmark-based label when MACRO_LABELS missing
label_lookup <- merge(label_lookup, top_hallmark[, .(macro_id = community_id, hallmark_term = term)],
                     by = "macro_id", all.x = TRUE)
label_lookup[is.na(label) | !(macro_id %in% names(MACRO_LABELS)),
             label := ifelse(!is.na(hallmark_term),
                             sub("HALLMARK_", "", hallmark_term),
                             paste0("Macro ", macro_id))]
label_lookup[macro_id == "Other", label := "Other"]

# ===========================================================================
# Panel A: Global Community Map
# ===========================================================================
# If a precomputed global_layout.csv exists, use it; otherwise fall back to
# computing a coarse layout from the composite edges (top edges only, to keep
# memory sane on 16K nodes x 3.6M edges).
# ===========================================================================
message("Computing / loading global layout ...")

layout_path <- file.path(NET_DIR, "global_layout.csv")

compute_layout_fallback <- function(top_edges, nodes_dt) {
  edges_dt <- top_edges[, .(gene_a, gene_b)]
  all_nodes <- unique(c(edges_dt$gene_a, edges_dt$gene_b))
  g <- graph_from_data_frame(edges_dt, directed = FALSE, vertices = data.frame(name = all_nodes))
  coords <- layout_with_fr(g, niter = 200, grid = "nogrid")
  data.table(gene = V(g)$name, x = coords[, 1], y = coords[, 2])
}

if (file.exists(layout_path)) {
  message("  Using precomputed ", layout_path)
  layout_dt <- fread(layout_path)
} else {
  message("  global_layout.csv missing -- computing FR layout from top composite edges")
  # Load only the strongest edges to keep layout tractable: top ~5% by p_composite
  comp_edges <- fread(file.path(NET_DIR, "composite_edges.csv"),
                      select = c("gene_a", "gene_b", "p_composite"))
  setorder(comp_edges, -p_composite)
  n_keep <- min(nrow(comp_edges), 5e5)
  top_edges <- comp_edges[1:n_keep]
  layout_dt <- compute_layout_fallback(top_edges, nodes)
  # Cache for future runs
  tryCatch(fwrite(layout_dt, layout_path),
           error = function(e) message("  (could not cache layout: ", conditionMessage(e), ")"))
  rm(comp_edges, top_edges); gc()
}

layout_dt <- merge(layout_dt, comm[, .(gene, macro_id, macro_label)], by = "gene", all.x = TRUE)
layout_dt <- merge(layout_dt, label_lookup, by.x = "macro_label", by.y = "macro_id", all.x = TRUE)
layout_dt[is.na(macro_label), macro_label := "Other"]
layout_dt[is.na(label), label := "Other"]

# Centroids for labels
centroids <- layout_dt[, .(cx = median(x, na.rm = TRUE),
                           cy = median(y, na.rm = TRUE),
                           n = .N),
                       by = .(macro_label, label)]
centroids <- centroids[n >= MIN_COMM_SIZE]

# Palette: fall back to hue_pal for any macro not in community_pal
all_macros <- unique(layout_dt$macro_label)
pal_a <- community_pal[names(community_pal) %in% all_macros]
missing_macros <- setdiff(all_macros, names(pal_a))
if (length(missing_macros) > 0) {
  extra <- scales::hue_pal()(length(missing_macros))
  names(extra) <- missing_macros
  pal_a <- c(pal_a, extra)
}

p_a <- ggplot(layout_dt, aes(x = x, y = y, color = macro_label)) +
  rasterize_layer(
    geom_point(size = 0.12, alpha = 0.35, shape = 16),
    dpi = 400
  ) +
  geom_mark_hull(aes(fill = macro_label, group = macro_label, filter = macro_label != "Other"),
                 concavity = 3, alpha = 0.07, linewidth = 0.2,
                 show.legend = FALSE) +
  geom_label_repel(data = centroids,
                   aes(x = cx, y = cy, label = label),
                   size = 2.0, fontface = "bold",
                   label.padding = 0.12, box.padding = 0.5,
                   segment.size = 0.15, max.overlaps = 30,
                   inherit.aes = FALSE, fill = "white", alpha = 0.9) +
  scale_color_manual(values = pal_a, guide = "none") +
  scale_fill_manual(values = pal_a, guide = "none") +
  labs(title = "Global community map (ForceAtlas2 layout)",
       subtitle = sprintf("%d genes in %d macro-communities (>=%d genes)",
                          nrow(layout_dt),
                          length(unique(layout_dt[macro_label != "Other", macro_label])),
                          MIN_COMM_SIZE),
       x = "FA2 dim 1", y = "FA2 dim 2") +
  theme_masld() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.subtitle = element_text(size = 6, color = "gray40"))

out_a <- file.path(OUT_DIR, "figS_network_final_panel_A_community_map.pdf")
save_fig(p_a, out_a, width = fig_full_width, height = 5.5)
message("Panel A saved: ", out_a)

# ===========================================================================
# Panel B: Modality Informativeness Heatmap
# ===========================================================================
# For each macro community, compute the mean per-layer posterior across edges
# where BOTH genes belong to the community. Use composite_edges.csv.
# ===========================================================================
message("Panel B: computing per-layer informativeness per community ...")

comp_edges <- fread(file.path(NET_DIR, "composite_edges.csv"))
# Attach community for both endpoints
gene_to_macro <- setNames(comm$macro_label, comm$gene)
comp_edges[, macro_a := gene_to_macro[gene_a]]
comp_edges[, macro_b := gene_to_macro[gene_b]]
within <- comp_edges[!is.na(macro_a) & !is.na(macro_b) & macro_a == macro_b]
within <- within[macro_a %in% keep_macros]

layer_cols <- paste0("p_", ACTIVE_LAYERS)
# Restrict to those actually present in composite_edges
layer_cols <- intersect(layer_cols, names(within))
active_layers_present <- sub("^p_", "", layer_cols)

# Mean posterior within each community x layer
if (nrow(within) > 0 && length(layer_cols) > 0) {
  info_mat <- within[, lapply(.SD, mean, na.rm = TRUE), by = macro_a,
                     .SDcols = layer_cols]
  mat_mat <- as.matrix(info_mat[, ..layer_cols])
  rownames(mat_mat) <- info_mat$macro_a
  colnames(mat_mat) <- active_layers_present

  # Order rows by size (largest community first)
  row_order <- comm_sizes[macro_id %in% rownames(mat_mat)][order(-N), macro_id]
  mat_mat <- mat_mat[row_order, , drop = FALSE]

  # Row labels = top hallmark / MACRO_LABELS
  row_labs <- sapply(rownames(mat_mat), function(m) {
    lab <- label_lookup[macro_id == m, label][1]
    if (is.na(lab)) m else paste0(m, ": ", lab)
  })

  col_labs <- layer_labels[colnames(mat_mat)]

  col_fun <- colorRamp2(
    c(0, quantile(mat_mat, 0.5, na.rm = TRUE), max(mat_mat, na.rm = TRUE)),
    c("#F7FBFF", "#6BAED6", "#08306B")
  )

  h <- Heatmap(mat_mat,
               name = "mean posterior",
               col = col_fun,
               cluster_rows = FALSE, cluster_columns = FALSE,
               row_labels = row_labs,
               column_labels = col_labs,
               row_names_side = "left",
               column_names_side = "bottom",
               column_names_rot = 45,
               row_names_gp = gpar(fontsize = 7),
               column_names_gp = gpar(fontsize = 7),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 v <- mat_mat[i, j]
                 if (!is.na(v) && v >= 0.01) {
                   grid.text(sprintf("%.2f", v), x, y,
                             gp = gpar(fontsize = 5,
                                       col = ifelse(v > 0.5, "white", "black")))
                 }
               },
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 6, fontface = "bold"),
                 labels_gp = gpar(fontsize = 5),
                 legend_height = unit(3, "cm")
               ),
               row_title = "Macro community",
               column_title = "Modality layer",
               row_title_gp = gpar(fontsize = 7, fontface = "bold"),
               column_title_gp = gpar(fontsize = 7, fontface = "bold"))

  out_b <- file.path(OUT_DIR, "figS_network_final_panel_B_informativeness.pdf")
  pdf(out_b, width = fig_col_width, height = 4.5)
  draw(h,
       heatmap_legend_side = "right",
       column_title = "Modality informativeness per community",
       column_title_gp = gpar(fontsize = 8, fontface = "bold"))
  dev.off()
  message("Panel B saved: ", out_b)
} else {
  message("Panel B: no within-community edges found; skipping")
}

rm(comp_edges, within); gc()

# ===========================================================================
# Panel C: Gold-Standard Recovery Benchmark
# ===========================================================================
message("Panel C: building gold-standard benchmark panel ...")

bench_path <- file.path(BENCH_DIR, "bench_per_layer.csv")
if (!file.exists(bench_path)) {
  message("  WARNING: ", bench_path, " not found; skipping Panel C")
} else {
  bench <- fread(bench_path)

  # Clean columns / coerce numerics
  for (col in c("fold_over_null", "fold_ci_lo", "fold_ci_hi", "p_perm")) {
    if (col %in% names(bench)) suppressWarnings(bench[, (col) := as.numeric(get(col))])
  }

  # Order layers: active layers first, composite variants last
  layer_order <- c(ACTIVE_LAYERS,
                   grep("composite", unique(bench$layer), value = TRUE))
  layer_order <- intersect(layer_order, unique(bench$layer))
  remaining <- setdiff(unique(bench$layer), layer_order)
  bench[, layer := factor(layer, levels = c(layer_order, remaining))]

  # Gold set labels
  gold_set_labels <- c(
    clinical_drug      = "Clinical drug targets",
    govaere25          = "Govaere 25-gene",
    resmetirom_pathway = "Resmetirom pathway"
  )
  bench[, gold_label := ifelse(gold_set %in% names(gold_set_labels),
                               gold_set_labels[gold_set],
                               gold_set)]

  # Cap CI for plotting (very wide CIs at small n dominate the x-axis)
  bench[, fold_plot    := pmin(fold_over_null, 50)]
  bench[, fold_ci_hi_p := pmin(fold_ci_hi, 50)]
  bench[, fold_ci_lo_p := pmax(fold_ci_lo, 0)]
  bench[, sig_star := fifelse(is.na(p_perm), "",
                     fifelse(p_perm < 0.001, "***",
                     fifelse(p_perm < 0.01,  "**",
                     fifelse(p_perm < 0.05,  "*", ""))))]

  # Composite highlight
  composite_layers <- grep("composite", levels(bench$layer), value = TRUE)
  bench[, is_composite := layer %in% composite_layers]

  p_c <- ggplot(bench, aes(x = layer, y = fold_plot, fill = is_composite)) +
    geom_col(width = 0.7, color = "gray20", linewidth = 0.2) +
    geom_errorbar(aes(ymin = fold_ci_lo_p, ymax = fold_ci_hi_p),
                  width = 0.25, linewidth = 0.25) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40",
               linewidth = 0.3) +
    geom_text(aes(label = sig_star, y = fold_ci_hi_p + 1.2),
              size = 2.2, vjust = 0) +
    facet_wrap(~ gold_label, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c("TRUE" = "#C2185B", "FALSE" = "#42A5F5"),
                      labels = c("TRUE" = "Composite", "FALSE" = "Single layer"),
                      name = NULL) +
    labs(title = "Gold-standard recovery: fold over null",
         subtitle = "Bars capped at 50x; dashed line = null (1x). * p<0.05, ** p<0.01, *** p<0.001",
         x = NULL, y = "Fold enrichment over null (+/- 95% CI)") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          legend.position = "top",
          plot.subtitle = element_text(size = 6, color = "gray40"),
          strip.text = element_text(size = 7, face = "bold"))

  out_c <- file.path(OUT_DIR, "figS_network_final_panel_C_benchmark.pdf")
  save_fig(p_c, out_c, width = fig_full_width, height = 3.5)
  message("Panel C saved: ", out_c)
}

# ===========================================================================
# Panel D: Example Gene Neighborhoods (THRB, HNF4A, PNPLA3, CFLAR)
# ===========================================================================
message("Panel D: building gene neighborhood panels ...")

example_genes <- c("THRB", "HNF4A", "PNPLA3", "CFLAR")
TOP_N_NEIGHBORS <- 15

build_neighborhood_plot <- function(gene_sym) {
  json_path <- file.path(GRAPH_DIR, paste0(gene_sym, ".json"))
  if (!file.exists(json_path)) {
    return(placeholder(paste0(gene_sym, "\n(JSON missing)")))
  }

  gd <- fromJSON(json_path, simplifyVector = FALSE)
  neigh <- gd$neighbors_composite
  if (is.null(neigh) || length(neigh) == 0) {
    return(placeholder(paste0(gene_sym, "\n(no neighbors)")))
  }

  # Flatten neighbor list into data.table: gene, p_composite, k_multiplicity, top_layer
  neigh_dt <- rbindlist(lapply(neigh, function(x) {
    layers_vec <- unlist(x$layers)
    if (length(layers_vec) == 0) {
      top_layer <- NA_character_
    } else {
      top_layer <- names(layers_vec)[which.max(layers_vec)]
    }
    data.table(
      gene          = x$gene,
      p_composite   = as.numeric(x$p_composite),
      k_multiplicity = as.integer(x$k_multiplicity),
      top_layer     = top_layer
    )
  }), fill = TRUE)

  setorder(neigh_dt, -p_composite)
  if (nrow(neigh_dt) > TOP_N_NEIGHBORS) neigh_dt <- neigh_dt[1:TOP_N_NEIGHBORS]

  # Build star graph: center -> each neighbor
  edges_df <- data.frame(
    from = gene_sym,
    to   = neigh_dt$gene,
    layer = neigh_dt$top_layer,
    p_composite = neigh_dt$p_composite,
    k_multiplicity = neigh_dt$k_multiplicity,
    stringsAsFactors = FALSE
  )
  nodes_df <- data.frame(
    name = c(gene_sym, neigh_dt$gene),
    is_center = c(TRUE, rep(FALSE, nrow(neigh_dt))),
    k_multiplicity = c(max(neigh_dt$k_multiplicity, na.rm = TRUE),
                       neigh_dt$k_multiplicity),
    stringsAsFactors = FALSE
  )
  g <- graph_from_data_frame(edges_df, vertices = nodes_df, directed = FALSE)

  # Community membership of center (for annotation)
  center_macro <- comm[gene == gene_sym, macro_id]
  center_label <- if (length(center_macro) == 0 || is.na(center_macro[1])) {
    "(unassigned)"
  } else {
    lab <- label_lookup[macro_id == center_macro[1], label][1]
    if (is.na(lab)) paste0("macro ", center_macro[1]) else lab
  }

  # Plot
  p_sub <- tryCatch({
    ggraph(g, layout = "fr") +
      geom_edge_link(aes(color = layer, width = p_composite),
                     alpha = 0.75) +
      geom_node_point(aes(size = k_multiplicity, fill = is_center),
                      shape = 21, stroke = 0.3, color = "gray30") +
      geom_node_text(aes(label = name,
                         fontface = ifelse(is_center, "bold.italic", "plain")),
                     size = 1.8, repel = TRUE, max.overlaps = 30,
                     bg.color = "white", bg.r = 0.08) +
      scale_edge_color_manual(values = layer_colors, na.value = "#BDBDBD",
                              name = "Top layer",
                              labels = layer_labels) +
      scale_edge_width_continuous(range = c(0.2, 1.1), guide = "none") +
      scale_fill_manual(values = c("TRUE" = "#FFD600", "FALSE" = "#ECEFF1"),
                        guide = "none") +
      scale_size_continuous(range = c(1.5, 4.5), name = "K_mult") +
      labs(title = sprintf("%s  -  %s", gene_sym, center_label)) +
      theme_masld() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.line = element_blank(), axis.title = element_blank(),
            legend.position = "none",
            plot.title = element_text(face = "bold", size = 7, hjust = 0.5))
  }, error = function(e) {
    message("  error plotting ", gene_sym, ": ", conditionMessage(e))
    placeholder(paste0(gene_sym, "\n(plot error)"))
  })

  p_sub
}

panels <- lapply(example_genes, build_neighborhood_plot)
names(panels) <- example_genes

# Shared legend: build once from a dummy data set that contains all layers
legend_df <- data.frame(
  x = 1, y = 1, layer = factor(names(layer_colors), levels = names(layer_colors))
)
legend_plot <- ggplot(legend_df, aes(x, y, color = layer)) +
  geom_point() +
  scale_color_manual(values = layer_colors, labels = layer_labels,
                     name = "Top modality layer") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        legend.text  = element_text(size = 6),
        legend.title = element_text(size = 7, face = "bold"))

extract_legend <- function(p) {
  gt <- ggplot_gtable(ggplot_build(p))
  idx <- which(sapply(gt$grobs, function(x) x$name == "guide-box"))
  if (length(idx) == 0) return(NULL)
  gt$grobs[[idx[1]]]
}
shared_leg <- extract_legend(legend_plot)

p_d_grid <- wrap_plots(panels, ncol = 2, nrow = 2)
p_d <- p_d_grid /
       wrap_elements(full = shared_leg) +
       plot_layout(heights = c(10, 1)) +
       plot_annotation(title = "Example gene neighborhoods (top-15 composite neighbors)",
                       theme = theme(plot.title = element_text(size = 8, face = "bold")))

out_d <- file.path(OUT_DIR, "figS_network_final_panel_D_neighborhoods.pdf")
save_fig(p_d, out_d, width = fig_full_width, height = 6.5)
message("Panel D saved: ", out_d)

elapsed <- (proc.time() - t0)["elapsed"]
message(sprintf("figS_network_final.R completed in %.1f seconds", elapsed))
