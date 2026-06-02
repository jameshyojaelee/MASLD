#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_271_prefig

# 271_preliminary_network_figures.R
# Quick preliminary figures from raw Stage 1 edges (pre-Bayesian)
# Shows: (A) edge inventory, (B) layer overlap, (C) example neighborhood, (D) degree distributions

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
})

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
NETDIR <- file.path(BASE, "RNA-seq/results/network")
OUTDIR <- file.path(BASE, "figures/misc/network_preliminary")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

LAYER_COLORS <- c(
  ppi = "#e74c3c", coexpr = "#3498db", regulon = "#2ecc71", lr = "#f39c12",
  genetic = "#9b59b6", pathway = "#1abc9c", spatial = "#e67e22", cosmos = "#34495e",
  cerna = "#e91e63", xspecies = "#607d8b"
)

LAYER_LABELS <- c(
  ppi = "PPI (STRING)", coexpr = "Co-expression", regulon = "TF Regulon",
  lr = "Ligand-Receptor", genetic = "Genetic (COLOC)", pathway = "Shared Pathway",
  spatial = "Spatial", cosmos = "COSMOS Path", cerna = "ceRNA", xspecies = "Cross-species"
)

nodes <- fread(file.path(NETDIR, "network_nodes.csv"))
cat("Nodes:", nrow(nodes), "\n")

layers <- c("ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "spatial", "cosmos", "cerna", "xspecies")
edge_list <- list()
for (l in layers) {
  f <- file.path(NETDIR, paste0("edges_", l, ".csv"))
  if (file.exists(f) && file.size(f) > 100) {
    dt <- fread(f)
    if (nrow(dt) > 0) {
      dt[, layer := l]
      edge_list[[l]] <- dt
      cat(l, ":", nrow(dt), "edges\n")
    }
  }
}

theme_net <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11)
  )

# ── Panel A: Edge Inventory (log-scale bar chart) ──
cat("\n--- Panel A: Edge Inventory ---\n")
inv <- data.table(
  layer = names(edge_list),
  n_edges = sapply(edge_list, nrow)
)
inv[, label := LAYER_LABELS[layer]]
inv[, label := factor(label, levels = label[order(n_edges)])]

pA <- ggplot(inv, aes(x = label, y = n_edges)) +
  geom_col(width = 0.7, fill = "#2c3e50") +
  geom_text(aes(label = formatC(n_edges, format = "d", big.mark = ",")),
            hjust = -0.1, size = 3.5) +
  scale_y_log10(labels = scales::comma, expand = expansion(mult = c(0, 0.3))) +
  coord_flip() +
  labs(
    title = "Multiplex Network: Edge Inventory",
    subtitle = paste0(nrow(nodes), " gene nodes across 10 evidence layers"),
    x = NULL, y = "Number of edges (log scale)"
  ) +
  theme_net

ggsave(file.path(OUTDIR, "panel_A_edge_inventory.pdf"), pA, width = 8, height = 5)
cat("Saved panel A\n")

# ── Panel B: Layer Overlap Heatmap ──
cat("\n--- Panel B: Layer Overlap ---\n")
gene_pairs <- list()
for (l in names(edge_list)) {
  dt <- edge_list[[l]]
  gene_pairs[[l]] <- paste(dt$gene_a, dt$gene_b, sep = "_")
}

active_layers <- names(gene_pairs)[sapply(gene_pairs, length) >= 10]
n_active <- length(active_layers)
jaccard_mat <- matrix(0, n_active, n_active, dimnames = list(active_layers, active_layers))
overlap_mat <- matrix(0, n_active, n_active, dimnames = list(active_layers, active_layers))

for (i in seq_along(active_layers)) {
  for (j in seq_along(active_layers)) {
    si <- gene_pairs[[active_layers[i]]]
    sj <- gene_pairs[[active_layers[j]]]
    inter <- length(intersect(si, sj))
    union_size <- length(union(si, sj))
    jaccard_mat[i, j] <- if (union_size > 0) inter / union_size else 0
    overlap_mat[i, j] <- inter
  }
}

jac_dt <- as.data.table(reshape2::melt(jaccard_mat, varnames = c("layer1", "layer2"), value.name = "jaccard"))
jac_dt[, label1 := LAYER_LABELS[as.character(layer1)]]
jac_dt[, label2 := LAYER_LABELS[as.character(layer2)]]

pB <- ggplot(jac_dt, aes(x = label1, y = label2, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(jaccard > 0.001, sprintf("%.3f", jaccard), "")),
            size = 3, color = "black") +
  scale_fill_gradient2(low = "white", mid = "#b3cde3", high = "#e74c3c", midpoint = 0.05,
                       name = "Jaccard\nIndex") +
  labs(
    title = "Edge Overlap Between Modality Layers",
    subtitle = "Jaccard index of gene-pair sets (low overlap = independent evidence)",
    x = NULL, y = NULL
  ) +
  theme_net +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTDIR, "panel_B_layer_overlap.pdf"), pB, width = 9, height = 7)
cat("Saved panel B\n")

# ── Panel C: Example Gene Neighborhoods (THRB, HNF4A, PNPLA3, CFLAR) ──
cat("\n--- Panel C: Example Neighborhoods ---\n")
example_genes <- c("THRB", "HNF4A", "PNPLA3", "CFLAR")

all_edges <- rbindlist(edge_list, use.names = TRUE, fill = TRUE)

plot_neighborhood <- function(gene, all_edges, top_k = 15) {
  neighbors <- all_edges[gene_a == gene | gene_b == gene]
  if (nrow(neighbors) == 0) return(NULL)

  neighbors[, partner := fifelse(gene_a == gene, gene_b, gene_a)]
  agg <- neighbors[, .(
    n_layers = uniqueN(layer),
    total_score = sum(raw_score, na.rm = TRUE),
    layers = paste(sort(unique(layer)), collapse = ",")
  ), by = partner]
  agg <- agg[order(-n_layers, -total_score)][1:min(top_k, .N)]

  sub_edges <- neighbors[partner %in% agg$partner]
  all_nodes <- unique(c(gene, agg$partner))

  g <- graph_from_data_frame(
    sub_edges[, .(from = fifelse(gene_a == gene, gene_a, gene_b),
                  to = fifelse(gene_a == gene, gene_b, gene_a),
                  layer = layer, weight = raw_score)],
    directed = FALSE, vertices = data.frame(name = all_nodes)
  )

  V(g)$is_center <- V(g)$name == gene
  V(g)$size <- ifelse(V(g)$is_center, 12, 6 + agg$n_layers[match(V(g)$name, agg$partner)])
  V(g)$color <- ifelse(V(g)$is_center, "#FFD700", "#CCCCCC")
  V(g)$label.cex <- ifelse(V(g)$is_center, 1.0, 0.7)
  E(g)$color <- LAYER_COLORS[E(g)$layer]

  return(g)
}

pdf(file.path(OUTDIR, "panel_C_example_neighborhoods.pdf"), width = 14, height = 12)
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
for (gene in example_genes) {
  g <- plot_neighborhood(gene, all_edges)
  if (!is.null(g)) {
    set.seed(42)
    plot(g, layout = layout_with_fr(g),
         vertex.label = V(g)$name,
         vertex.size = V(g)$size,
         vertex.color = V(g)$color,
         vertex.frame.color = NA,
         edge.color = adjustcolor(E(g)$color, alpha.f = 0.7),
         edge.width = 1.5,
         main = paste0(gene, " neighborhood (top 15)"))
  } else {
    plot.new()
    title(paste0(gene, " — no edges found"))
  }
}

plot.new()
legend("center", legend = LAYER_LABELS[names(LAYER_COLORS)],
       fill = LAYER_COLORS, border = NA, cex = 0.9, ncol = 2,
       title = "Edge Layer")
dev.off()
cat("Saved panel C\n")

# ── Panel D: Degree Distribution Per Layer ──
cat("\n--- Panel D: Degree Distributions ---\n")
degree_list <- list()
for (l in names(edge_list)) {
  dt <- edge_list[[l]]
  genes <- c(dt$gene_a, dt$gene_b)
  deg <- as.data.table(table(genes))
  setnames(deg, c("gene", "degree"))
  deg[, layer := l]
  degree_list[[l]] <- deg
}
all_deg <- rbindlist(degree_list)
all_deg[, label := LAYER_LABELS[layer]]

pD <- ggplot(all_deg[layer %in% c("ppi", "coexpr", "pathway", "lr", "genetic")],
       aes(x = degree, fill = layer, color = layer)) +
  geom_density(alpha = 0.3, adjust = 2) +
  scale_x_log10() +
  scale_fill_manual(values = LAYER_COLORS, labels = LAYER_LABELS) +
  scale_color_manual(values = LAYER_COLORS, labels = LAYER_LABELS) +
  labs(
    title = "Degree Distribution by Layer",
    subtitle = "Most genes have few connections; hub genes emerge in PPI and pathway layers",
    x = "Degree (log scale)", y = "Density",
    fill = "Layer", color = "Layer"
  ) +
  theme_net

ggsave(file.path(OUTDIR, "panel_D_degree_distributions.pdf"), pD, width = 9, height = 5)
cat("Saved panel D\n")

# ── Panel E: Multi-layer Gene Connectivity ──
cat("\n--- Panel E: Multi-layer connectivity ---\n")
gene_layer_counts <- all_edges[, .(n_layers = uniqueN(layer)), by = .(gene_a, gene_b)]
pair_layer_dist <- gene_layer_counts[, .N, by = n_layers][order(n_layers)]
pair_layer_dist[, pct := N / sum(N) * 100]

pE <- ggplot(pair_layer_dist, aes(x = factor(n_layers), y = N)) +
  geom_col(fill = "#2c3e50", width = 0.5) +
  geom_text(aes(label = paste0(formatC(N, big.mark = ","), "\n(", round(pct, 1), "%)")),
            vjust = -0.2, size = 3.5) +
  scale_y_log10(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Gene Pair Connectivity Across Layers",
    subtitle = "How many modalities independently connect the same gene pair?",
    x = "Number of layers connecting the pair", y = "Number of gene pairs (log scale)"
  ) +
  theme_net

ggsave(file.path(OUTDIR, "panel_E_multilayer_connectivity.pdf"), pE, width = 4.5, height = 5)
cat("Saved panel E\n")

# ── Panel F: Hub Genes Table ──
cat("\n--- Panel F: Hub gene summary ---\n")
gene_degree <- all_edges[, .(
  total_edges = .N,
  n_layers = uniqueN(layer),
  layers = paste(sort(unique(layer)), collapse = ", ")
), by = .(gene = gene_a)]

gene_degree2 <- all_edges[, .(
  total_edges = .N,
  n_layers = uniqueN(layer),
  layers = paste(sort(unique(layer)), collapse = ", ")
), by = .(gene = gene_b)]

combined_deg <- rbind(gene_degree, gene_degree2)[,
  .(total_edges = sum(total_edges), n_layers = max(n_layers)),
  by = gene][order(-n_layers, -total_edges)]

hub_genes <- merge(combined_deg[1:30], nodes[, .(human_symbol, is_deg, is_conserved,
  coloc_susie_best_pp4, sex_class, dgidb_druggable)],
  by.x = "gene", by.y = "human_symbol", all.x = TRUE)

fwrite(hub_genes, file.path(OUTDIR, "hub_genes_top30.csv"))
cat("Top 30 hub genes saved\n")
cat("\nTop 10 hub genes:\n")
print(hub_genes[1:10, .(gene, total_edges, n_layers, is_deg, coloc = round(coloc_susie_best_pp4, 3), druggable = dgidb_druggable)])

cat("\n\nElapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "mins\n")
cat("Output:", OUTDIR, "\n")
cat("Done.\n")
