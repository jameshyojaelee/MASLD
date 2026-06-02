#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_272_bench

# 272_network_benchmark.R
# Benchmark multi-evidence composite network vs single-layer baselines
# Response to adversarial review (Tier 3 fix):
#   - Does composite recover MASLD gold-standard gene sets better than
#     any single layer or a random network?
#
# Gold standards:
#   1. Positive control genes (65 curated MASLD-relevant genes)
#   2. Clinical drug target genes (from positive_control.csv)
#   3. DEGs + COLOC genes (multi-evidence gold standard)
#
# Metrics:
#   - Enrichment of gold-standard genes among high-degree nodes per layer
#   - Connectivity within gold-standard set (fraction of pairs with edges)
#   - Precision-recall for gold-standard link prediction

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
})

t0 <- Sys.time()

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
NETDIR <- file.path(BASE, "RNA-seq/results/network")
OUTDIR <- file.path(BASE, "figures/misc/network_benchmark")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "data"), showWarnings = FALSE)

nodes <- fread(file.path(NETDIR, "network_nodes.csv"))
cat("Nodes:", nrow(nodes), "\n")

# --- Gold standards ---
pc <- fread(file.path(BASE, "results/library/positive_control.csv"))
cat("Positive controls:", nrow(pc), "\n")

gold_pos_ctrl <- intersect(pc$`Gene symbol`, nodes$human_symbol)
cat("  Positive controls in V:", length(gold_pos_ctrl), "\n")

gold_clinical <- pc[!is.na(Clinical_Drug) & Clinical_Drug != "", `Gene symbol`]
gold_clinical <- intersect(gold_clinical, nodes$human_symbol)
cat("  Clinical drug targets in V:", length(gold_clinical), "\n")

gold_coloc <- nodes[!is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.5, human_symbol]
cat("  COLOC genes in V:", length(gold_coloc), "\n")

gold_deg_coloc <- nodes[!is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.5 & is_deg == TRUE, human_symbol]
cat("  DEG + COLOC genes in V:", length(gold_deg_coloc), "\n")

gold_sets <- list(
  positive_control = gold_pos_ctrl,
  clinical_drug    = gold_clinical,
  coloc            = gold_coloc,
  deg_and_coloc    = gold_deg_coloc
)

# --- Load layer edges ---
layers <- c("ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna")
edge_files <- list()
for (l in layers) {
  f <- file.path(NETDIR, paste0("edges_", l, ".csv"))
  if (!file.exists(f)) next
  dt <- fread(f, select = c("gene_a", "gene_b", "raw_score"))
  dt[, layer := l]
  edge_files[[l]] <- dt
  cat(l, ":", nrow(dt), "edges\n")
}

# Composite at K >= 2 (the key test: does convergence tell us more than any single layer?)
comp_path <- file.path(NETDIR, "composite_edges.csv")
if (file.exists(comp_path)) {
  comp <- fread(comp_path, select = c("gene_a", "gene_b", "p_composite", "k_multiplicity"))
  comp_k2 <- comp[k_multiplicity >= 2]
  comp_k3 <- comp[k_multiplicity >= 3]
  edge_files[["composite_k2"]] <- comp_k2[, .(gene_a, gene_b, raw_score = p_composite, layer = "composite_k2")]
  edge_files[["composite_k3"]] <- comp_k3[, .(gene_a, gene_b, raw_score = p_composite, layer = "composite_k3")]
  edge_files[["composite_p09"]] <- comp[p_composite > 0.9,
                                        .(gene_a, gene_b, raw_score = p_composite, layer = "composite_p09")]
  cat("composite_k2 :", nrow(comp_k2), "edges\n")
  cat("composite_k3 :", nrow(comp_k3), "edges\n")
}

# --- Benchmark 1: Gold-standard connectivity per layer ---
# For each layer, what fraction of gold-standard gene pairs are connected?

compute_gs_connectivity <- function(el, gold_genes) {
  if (length(gold_genes) < 2) return(NA)
  total_pairs <- choose(length(gold_genes), 2)
  gs <- el[gene_a %in% gold_genes & gene_b %in% gold_genes]
  nrow(gs) / total_pairs
}

bench1 <- data.table()
for (ln in names(edge_files)) {
  for (gn in names(gold_sets)) {
    frac <- compute_gs_connectivity(edge_files[[ln]], gold_sets[[gn]])
    bench1 <- rbind(bench1, data.table(
      layer = ln, gold_set = gn, n_gold = length(gold_sets[[gn]]),
      fraction_connected = frac
    ))
  }
}

fwrite(bench1, file.path(OUTDIR, "data", "bench1_gs_connectivity.csv"))

# --- Benchmark 2: Degree enrichment for gold-standard genes ---
# Are gold-standard genes HIGHER degree than random nodes?

compute_degree_enrichment <- function(el, gold_genes, all_genes) {
  deg_dt <- rbind(
    el[, .N, by = gene_a][, .(gene = gene_a, deg = N)],
    el[, .N, by = gene_b][, .(gene = gene_b, deg = N)]
  )[, .(deg = sum(deg)), by = gene]

  gs_deg <- deg_dt[gene %in% gold_genes, deg]
  bg_deg <- deg_dt[!(gene %in% gold_genes), deg]

  if (length(gs_deg) < 5 || length(bg_deg) < 5) return(c(NA, NA, NA))
  w <- wilcox.test(gs_deg, bg_deg, alternative = "greater")
  c(mean_gs = mean(gs_deg), mean_bg = mean(bg_deg), p_wilcox = w$p.value)
}

bench2 <- data.table()
for (ln in names(edge_files)) {
  for (gn in names(gold_sets)) {
    r <- compute_degree_enrichment(edge_files[[ln]], gold_sets[[gn]], nodes$human_symbol)
    bench2 <- rbind(bench2, data.table(
      layer = ln, gold_set = gn,
      mean_deg_gold = r[1], mean_deg_background = r[2], p_wilcox = r[3]
    ))
  }
}

fwrite(bench2, file.path(OUTDIR, "data", "bench2_degree_enrichment.csv"))

# --- Benchmark 3: Random baseline ---
# Compare to N random edge sets of the same size per layer
set.seed(42)
n_perm <- 100

random_baseline <- function(n_edges, gold_genes, all_genes) {
  n_v <- length(all_genes)
  fracs <- sapply(seq_len(n_perm), function(i) {
    ga <- sample(all_genes, n_edges, replace = TRUE)
    gb <- sample(all_genes, n_edges, replace = TRUE)
    sub <- data.table(gene_a = ga, gene_b = gb)
    sub <- sub[gene_a != gene_b]
    total_pairs <- choose(length(gold_genes), 2)
    n_gs <- nrow(sub[gene_a %in% gold_genes & gene_b %in% gold_genes])
    n_gs / total_pairs
  })
  c(mean = mean(fracs), sd = sd(fracs), p95 = quantile(fracs, 0.95))
}

bench3 <- data.table()
for (ln in names(edge_files)) {
  n_edges <- nrow(edge_files[[ln]])
  for (gn in names(gold_sets)) {
    observed <- compute_gs_connectivity(edge_files[[ln]], gold_sets[[gn]])
    rand <- random_baseline(n_edges, gold_sets[[gn]], nodes$human_symbol)
    enrichment <- if (!is.na(observed) && !is.na(rand[1]) && rand[1] > 0) observed / rand[1] else NA
    bench3 <- rbind(bench3, data.table(
      layer = ln, gold_set = gn,
      observed = observed, random_mean = rand[1], random_sd = rand[2],
      random_p95 = rand[3], enrichment_fold = enrichment
    ))
  }
}

fwrite(bench3, file.path(OUTDIR, "data", "bench3_random_baseline.csv"))

# --- Figure: enrichment fold per layer × gold set ---
bench3_plot <- bench3[!is.na(enrichment_fold) & enrichment_fold > 0]
bench3_plot[, layer := factor(layer, levels = c(
  "ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna",
  "composite_p09", "composite_k2", "composite_k3"
))]
bench3_plot[, gold_set := factor(gold_set, levels = c(
  "positive_control", "clinical_drug", "coloc", "deg_and_coloc"
))]

p_fold <- ggplot(bench3_plot, aes(x = gold_set, y = layer, fill = log10(enrichment_fold))) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1fx", enrichment_fold)), size = 3.5, color = "black") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, name = "log10\nfold-enrichment\nvs random") +
  labs(
    title = "Gold-standard recovery: per-layer fold-enrichment over random baseline",
    subtitle = "Higher = the layer connects gold-standard genes more than expected by chance",
    x = "Gold-standard gene set", y = "Layer / composite"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(OUTDIR, "bench_fold_enrichment.pdf"), p_fold, width = 8, height = 5)

# --- Summary: does composite beat single layers? ---
summary_row <- function(layer_name) {
  bench3[layer == layer_name, .(mean_enrichment = mean(enrichment_fold, na.rm = TRUE))]
}

composite_beats <- data.table(
  layer = unique(bench3$layer),
  mean_fold = sapply(unique(bench3$layer), function(l) mean(bench3[layer == l, enrichment_fold], na.rm = TRUE))
)[order(-mean_fold)]

fwrite(composite_beats, file.path(OUTDIR, "data", "layer_ranking.csv"))
cat("\n=== Layer ranking by mean fold-enrichment ===\n")
print(composite_beats)

cat("\nElapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "mins\n")
cat("Outputs:", OUTDIR, "\n")
