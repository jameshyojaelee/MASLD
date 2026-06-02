#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --job-name=net_272b_bench
#SBATCH --output=logs/net_272b_bench_%j.out
#SBATCH --error=logs/net_272b_bench_%j.err
#
# 272b_robust_benchmark.R
# Rigorous re-benchmark replacing the flawed 272 analysis.
#
# Fixes over 272:
#   (1) Degree-preserving Maslov-Sneppen null (igraph::rewire / keeping_degseq)
#       instead of IID random-edge sampling that inflates enrichment.
#   (2) 1,000 rewirings per layer for tight 95% CIs (was 100).
#   (3) Non-tautological, (mostly) held-out gold standards. Tautological sets
#       like "COLOC genes in the genetic layer" are REMOVED.
#   (4) ROC/PR: composite vs. best single layer for gold-pair recovery.
#
# Environment: micromamba activate rnaseq
# Depends on: 273_normalize_ids.py must have completed (produces
#             results/network/id_normalization_summary.csv)
#
# Outputs:
#   figures/misc/network_benchmark_v2/data/bench_per_layer.csv
#   figures/misc/network_benchmark_v2/bench_fold_enrichment_heatmap.pdf
#   figures/misc/network_benchmark_v2/bench_ROC_composite_vs_layers.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(pROC)
  library(parallel)
})

t0 <- Sys.time()

# ---------------------------------------------------------------------------
# Paths & config
# ---------------------------------------------------------------------------
BASE   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
NETDIR <- file.path(BASE, "RNA-seq/results/network")
OUTDIR <- file.path(BASE, "figures/misc/network_benchmark_v2")
DATDIR <- file.path(OUTDIR, "data")
dir.create(DATDIR, recursive = TRUE, showWarnings = FALSE)

N_REWIRE    <- 1000    # null rewirings per layer
SWAP_MULT   <- 10      # edge swaps per rewiring = SWAP_MULT * n_edges
N_CORES     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
if (is.na(N_CORES) || N_CORES < 1) N_CORES <- 8
MC_CORES    <- min(8L, N_CORES)   # rewiring is single-threaded igraph; 8 parallel tasks
SEED        <- 42
set.seed(SEED)

cat("======================================================================\n")
cat("272b ROBUST NETWORK BENCHMARK\n")
cat("  Cores available: ", N_CORES, "  | mclapply cores: ", MC_CORES, "\n", sep = "")
cat("  N_REWIRE: ", N_REWIRE, "  | SWAP_MULT: ", SWAP_MULT, "\n", sep = "")
cat("======================================================================\n\n")

# ---------------------------------------------------------------------------
# Guard: require 273 to have completed
# ---------------------------------------------------------------------------
id_norm_path <- file.path(NETDIR, "id_normalization_summary.csv")
if (!file.exists(id_norm_path)) {
  stop("Run 273_normalize_ids.py first (missing: ", id_norm_path, ")")
}
cat("ID normalization summary present: ", id_norm_path, "\n\n", sep = "")

# ---------------------------------------------------------------------------
# Load nodes
# ---------------------------------------------------------------------------
nodes <- fread(file.path(NETDIR, "network_nodes.csv"))
all_genes <- unique(nodes$human_symbol)
cat("Nodes in V: ", length(all_genes), "\n", sep = "")

# ---------------------------------------------------------------------------
# Gold standards (non-tautological)
# ---------------------------------------------------------------------------
pc_path <- file.path(BASE, "results/library/positive_control.csv")
pc <- fread(pc_path)

# (a) clinical_drug: drugs-approved MASLD targets. Orthogonal to PPI / coexpr /
#     pathway / regulon / cerna. Could overlap with genetic (if drug target is
#     also a GWAS locus) -- document this, do not hide it.
clin_col <- intersect(c("Clinical_Drug", "clinical_drug"), names(pc))
if (length(clin_col) > 0) {
  gold_clinical <- pc[!is.na(get(clin_col[1])) & get(clin_col[1]) != "", `Gene symbol`]
} else {
  gold_clinical <- character(0)
}
gold_clinical <- intersect(gold_clinical, all_genes)

# (b) Govaere 25-gene prognostic panel (Govaere et al. 2020, Sci Transl Med).
#     Orthogonal to ALL construction sources (panel defined by bulk RNA-seq
#     prognostic modelling; none of our layers "know" this panel).
govaere25_published <- c(
  "AKR1B10", "SOX4", "TTC39A", "VIL1", "CLIC6", "ABCB4", "LTBP2", "STMN2",
  "AGMO", "ANGPTL4", "CA3", "CCBE1", "CFB", "CSAD", "FABP4", "FADS1",
  "IGFBP7", "MME", "PANX2", "PLS1", "RP11-1399P15.1", "SDC1", "SRPX",
  "TMEM9", "TPBG"
)
val_dir <- file.path(BASE, "RNA-seq/results/validation")
govaere25 <- govaere25_published
# Try to load a richer file if present
for (fname in c("govaere25.csv", "govaere_panel.csv", "published_panels.csv")) {
  fpath <- file.path(val_dir, fname)
  if (file.exists(fpath)) {
    try({
      tmp <- fread(fpath)
      sym_col <- intersect(c("Gene symbol", "gene_symbol", "symbol", "gene"),
                           names(tmp))
      if (length(sym_col) > 0) {
        govaere25 <- unique(c(govaere25, tmp[[sym_col[1]]]))
      }
    }, silent = TRUE)
    break
  }
}
gold_govaere <- intersect(govaere25, all_genes)

# (c) resmetirom_pathway: THRB + canonical downstream targets. Orthogonal to
#     all statistical layers (it's a pharmacological/mechanistic set).
resm_pathway <- c("THRB", "PCSK9", "PPARA", "APOA1")
gold_resm <- intersect(resm_pathway, all_genes)

# (d) opentargets_masld (optional)
ot_path <- file.path(BASE, "data/opentargets_masld.csv")
gold_ot <- character(0)
if (file.exists(ot_path)) {
  ot <- fread(ot_path)
  sym_col <- intersect(c("Gene symbol", "gene_symbol", "symbol", "gene"), names(ot))
  if (length(sym_col) > 0) {
    gold_ot <- intersect(unique(ot[[sym_col[1]]]), all_genes)
  }
} else {
  cat("WARNING: OpenTargets gold standard missing at ", ot_path, " -- skipping\n",
      sep = "")
}

gold_sets <- list(
  clinical_drug      = gold_clinical,
  govaere25          = gold_govaere,
  resmetirom_pathway = gold_resm
)
if (length(gold_ot) >= 2) gold_sets[["opentargets_masld"]] <- gold_ot

# Print gold-set summary and overlaps
cat("\nGold-standard sizes (in V):\n")
for (gn in names(gold_sets)) {
  cat(sprintf("  %-22s : %4d\n", gn, length(gold_sets[[gn]])))
}
cat("\nPairwise overlap between gold sets:\n")
gn_all <- names(gold_sets)
for (i in seq_along(gn_all)) {
  for (j in seq_along(gn_all)) {
    if (i >= j) next
    ov <- length(intersect(gold_sets[[gn_all[i]]], gold_sets[[gn_all[j]]]))
    cat(sprintf("  %-22s ∩ %-22s = %d\n", gn_all[i], gn_all[j], ov))
  }
}
cat("\n")

# ---------------------------------------------------------------------------
# Load layer edges
# ---------------------------------------------------------------------------
layers <- c("ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna")
edge_files <- list()
for (l in layers) {
  f <- file.path(NETDIR, paste0("edges_", l, ".csv"))
  if (!file.exists(f)) {
    cat("  (missing: ", l, ")\n", sep = "")
    next
  }
  dt <- fread(f, select = c("gene_a", "gene_b"))
  # Deduplicate undirected edges (canonical a<b ordering)
  dt[, c("gene_a", "gene_b") := list(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
  dt <- unique(dt[gene_a != gene_b & gene_a %in% all_genes & gene_b %in% all_genes])
  edge_files[[l]] <- dt
  cat(sprintf("  %-14s : %7d edges (deduped, in V)\n", l, nrow(dt)))
}

# Composite slices
comp_path <- file.path(NETDIR, "composite_edges.csv")
if (file.exists(comp_path)) {
  comp <- fread(comp_path, select = c("gene_a", "gene_b", "p_composite", "k_multiplicity"))
  comp[, c("gene_a", "gene_b") := list(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
  comp <- unique(comp[gene_a != gene_b & gene_a %in% all_genes & gene_b %in% all_genes])

  edge_files[["composite_p05"]] <- comp[p_composite > 0.5,  .(gene_a, gene_b)]
  edge_files[["composite_p07"]] <- comp[p_composite > 0.7,  .(gene_a, gene_b)]
  edge_files[["composite_p09"]] <- comp[p_composite > 0.9,  .(gene_a, gene_b)]
  edge_files[["composite_p095"]] <- comp[p_composite > 0.95, .(gene_a, gene_b)]
  edge_files[["composite_k2"]]  <- comp[k_multiplicity >= 2, .(gene_a, gene_b)]
  edge_files[["composite_k3"]]  <- comp[k_multiplicity >= 3, .(gene_a, gene_b)]

  for (cn in c("composite_p05", "composite_p07", "composite_p09", "composite_p095",
               "composite_k2", "composite_k3")) {
    cat(sprintf("  %-14s : %7d edges\n", cn, nrow(edge_files[[cn]])))
  }
}
cat("\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Observed fraction of gold-pairs connected by edges in `el`
gs_connectivity <- function(el, gold_genes) {
  if (length(gold_genes) < 2 || nrow(el) == 0) return(NA_real_)
  n_pairs_total <- choose(length(gold_genes), 2)
  gs <- el[gene_a %in% gold_genes & gene_b %in% gold_genes]
  nrow(gs) / n_pairs_total
}

# Degree enrichment: Wilcoxon U-stat comparing degree of gold genes vs non-gold
degree_enrichment <- function(el, gold_genes) {
  if (nrow(el) == 0 || length(gold_genes) < 5) return(list(U = NA, p = NA, mean_g = NA, mean_bg = NA))
  deg_dt <- rbind(
    el[, .(gene = gene_a)],
    el[, .(gene = gene_b)]
  )[, .(deg = .N), by = gene]
  # Add zero-degree genes
  missing_g <- setdiff(all_genes, deg_dt$gene)
  if (length(missing_g) > 0) {
    deg_dt <- rbind(deg_dt, data.table(gene = missing_g, deg = 0L))
  }
  gs_deg <- deg_dt[gene %in% gold_genes, deg]
  bg_deg <- deg_dt[!(gene %in% gold_genes), deg]
  if (length(gs_deg) < 5 || length(bg_deg) < 5)
    return(list(U = NA, p = NA, mean_g = NA, mean_bg = NA))
  w <- suppressWarnings(wilcox.test(gs_deg, bg_deg, alternative = "greater"))
  list(U = unname(w$statistic), p = w$p.value,
       mean_g = mean(gs_deg), mean_bg = mean(bg_deg))
}

# Build igraph from edge list (undirected, simple)
build_graph <- function(el) {
  if (nrow(el) == 0) return(NULL)
  g <- graph_from_data_frame(el[, .(gene_a, gene_b)], directed = FALSE,
                             vertices = data.frame(name = all_genes))
  simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

# Maslov-Sneppen degree-preserving rewiring via igraph (C impl)
rewire_one <- function(g, n_swaps, gold_genes, seed_i) {
  set.seed(seed_i)
  gr <- rewire(g, keeping_degseq(loops = FALSE, niter = n_swaps))
  ed <- as_edgelist(gr, names = TRUE)
  if (nrow(ed) == 0) return(NA_real_)
  a <- pmin(ed[, 1], ed[, 2]); b <- pmax(ed[, 1], ed[, 2])
  gs_mask <- a %in% gold_genes & b %in% gold_genes
  n_pairs_total <- choose(length(gold_genes), 2)
  sum(gs_mask) / n_pairs_total
}

# ---------------------------------------------------------------------------
# Main loop: per (layer, gold_set), compute observed + null distribution
# ---------------------------------------------------------------------------
layer_names <- names(edge_files)
results <- list()

for (ln in layer_names) {
  el <- edge_files[[ln]]
  t_layer <- Sys.time()
  cat("==== Layer: ", ln, " (", nrow(el), " edges) ====\n", sep = "")

  if (nrow(el) == 0) {
    cat("  (skipping: 0 edges)\n")
    next
  }

  g <- build_graph(el)
  n_edges <- ecount(g)
  n_swaps <- SWAP_MULT * n_edges

  for (gn in names(gold_sets)) {
    gold <- gold_sets[[gn]]
    if (length(gold) < 2) {
      cat("  ", gn, ": < 2 gold genes, skipping\n", sep = "")
      next
    }

    observed <- gs_connectivity(el, gold)
    de       <- degree_enrichment(el, gold)

    # Null: N_REWIRE degree-preserving rewirings, parallel
    t_null <- Sys.time()
    null_fracs <- unlist(mclapply(seq_len(N_REWIRE), function(i) {
      rewire_one(g, n_swaps, gold, seed_i = SEED + 1000L * which(layer_names == ln) + i)
    }, mc.cores = MC_CORES, mc.preschedule = TRUE))
    null_time <- as.numeric(difftime(Sys.time(), t_null, units = "secs"))

    mean_null <- mean(null_fracs, na.rm = TRUE)
    sd_null   <- sd(null_fracs,   na.rm = TRUE)
    fold      <- if (!is.na(mean_null) && mean_null > 0) observed / mean_null else NA_real_

    # p-value: fraction of null >= observed (add +1/+1 smoothing)
    n_ge <- sum(null_fracs >= observed, na.rm = TRUE)
    p_perm <- (n_ge + 1) / (sum(!is.na(null_fracs)) + 1)

    # 95% CI on fold: use null bootstrap -> observed / resampled_null_mean
    # Here we report the 95% range of observed / null_i as the fold CI.
    fold_dist <- observed / null_fracs
    fold_dist <- fold_dist[is.finite(fold_dist)]
    ci_lo <- if (length(fold_dist) > 10) quantile(fold_dist, 0.025, na.rm = TRUE) else NA_real_
    ci_hi <- if (length(fold_dist) > 10) quantile(fold_dist, 0.975, na.rm = TRUE) else NA_real_

    results[[length(results) + 1L]] <- data.table(
      layer             = ln,
      gold_set          = gn,
      n_gold            = length(gold),
      n_edges           = n_edges,
      fraction_connected = observed,
      null_mean         = mean_null,
      null_sd           = sd_null,
      fold_over_null    = fold,
      fold_ci_lo        = ci_lo,
      fold_ci_hi        = ci_hi,
      p_perm            = p_perm,
      n_null            = sum(!is.na(null_fracs)),
      degree_U          = de$U,
      degree_p          = de$p,
      mean_deg_gold     = de$mean_g,
      mean_deg_bg       = de$mean_bg,
      null_seconds      = null_time
    )

    cat(sprintf("  %-22s obs=%.4g  null=%.4g  fold=%.2fx [%.2f-%.2f]  p=%.1e  (%.0fs)\n",
                gn, observed, mean_null, fold,
                ifelse(is.na(ci_lo), NA_real_, ci_lo),
                ifelse(is.na(ci_hi), NA_real_, ci_hi),
                p_perm, null_time))
  }
  cat("  [layer elapsed: ",
      round(as.numeric(difftime(Sys.time(), t_layer, units = "mins")), 2),
      " min]\n\n", sep = "")
}

bench <- rbindlist(results)
fwrite(bench, file.path(DATDIR, "bench_per_layer.csv"))
cat("\nWrote: ", file.path(DATDIR, "bench_per_layer.csv"), "\n\n", sep = "")

# ---------------------------------------------------------------------------
# Figure 1: fold-enrichment heatmap with 95% CI text
# ---------------------------------------------------------------------------
bench_plot <- copy(bench)
bench_plot <- bench_plot[!is.na(fold_over_null) & is.finite(fold_over_null) &
                           fold_over_null > 0]

layer_order <- c(
  "ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna",
  "composite_p05", "composite_p07", "composite_p09", "composite_p095",
  "composite_k2", "composite_k3"
)
layer_order <- intersect(layer_order, bench_plot$layer)
bench_plot[, layer := factor(layer, levels = layer_order)]

gs_order <- intersect(
  c("clinical_drug", "govaere25", "resmetirom_pathway", "opentargets_masld"),
  unique(bench_plot$gold_set)
)
bench_plot[, gold_set := factor(gold_set, levels = gs_order)]

bench_plot[, label := sprintf("%.1fx\n[%.1f-%.1f]",
                              fold_over_null, fold_ci_lo, fold_ci_hi)]
bench_plot[!is.finite(fold_ci_lo) | !is.finite(fold_ci_hi),
           label := sprintf("%.1fx", fold_over_null)]
bench_plot[, log10_fold := log10(pmax(fold_over_null, 1e-3))]

p_heat <- ggplot(bench_plot, aes(x = gold_set, y = layer, fill = log10_fold)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 2.8, color = "black", lineheight = 0.9) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, name = "log10 fold\n(deg-pres null)") +
  labs(
    title = "Gold-standard recovery: fold-enrichment vs. degree-preserving null",
    subtitle = sprintf("N=%d Maslov-Sneppen rewirings per layer; 95%% CI in brackets",
                       N_REWIRE),
    x = "Gold-standard gene set", y = "Layer / composite"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid  = element_blank())

ggsave(file.path(OUTDIR, "bench_fold_enrichment_heatmap.pdf"),
       p_heat, width = 9, height = 6)
cat("Wrote: ", file.path(OUTDIR, "bench_fold_enrichment_heatmap.pdf"), "\n", sep = "")

# ---------------------------------------------------------------------------
# Figure 2: ROC / PR -- composite vs best single layer
# ---------------------------------------------------------------------------
# Task: for each gold set, for each candidate "scoring scheme" (layer), rank
# all gene-pairs (restricted to gene_a, gene_b in V) by edge presence, and ask
# how well it recovers gold-pair connections.
#
# Positives = gold-gold pairs (both endpoints in gold set)
# Negatives = gold-nongold or nongold-nongold pairs (sampled for tractability)
#
# Because the positive set is small, we can enumerate all gold-gold pairs
# and sample a matched-size negative set of gold-nongold pairs (realistic
# "distractors"), then score each pair per layer as (edge present ? w : 0).

roc_layers_single <- c("ppi", "coexpr", "regulon", "lr", "genetic", "pathway", "cerna")
roc_layers_single <- intersect(roc_layers_single, names(edge_files))

# Composite scoring: use p_composite directly if composite loaded
has_composite <- file.exists(comp_path)
if (has_composite) {
  comp_dt <- fread(comp_path, select = c("gene_a", "gene_b", "p_composite"))
  comp_dt[, c("gene_a", "gene_b") := list(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
  comp_dt <- unique(comp_dt[gene_a != gene_b & gene_a %in% all_genes & gene_b %in% all_genes])
  setkey(comp_dt, gene_a, gene_b)
}

# Edge lookup per layer as keyed data.tables
layer_dt <- lapply(roc_layers_single, function(ln) {
  dt <- edge_files[[ln]][, .(gene_a, gene_b)]
  dt[, in_layer := 1L]
  setkey(dt, gene_a, gene_b)
  dt
})
names(layer_dt) <- roc_layers_single

roc_rows <- list()
pr_rows  <- list()

for (gn in names(gold_sets)) {
  gold <- gold_sets[[gn]]
  if (length(gold) < 3) next

  # Enumerate gold-gold pairs (positives)
  gg <- as.data.table(t(combn(sort(gold), 2)))
  setnames(gg, c("gene_a", "gene_b"))
  gg[, label := 1L]

  # Sample gold-nongold pairs as negatives (~5x positives, capped)
  n_pos <- nrow(gg)
  n_neg <- min(50000L, max(5L * n_pos, 5000L))
  nongold <- setdiff(all_genes, gold)
  set.seed(SEED + which(names(gold_sets) == gn))
  ng <- data.table(
    gene_a = sample(gold,    n_neg, replace = TRUE),
    gene_b = sample(nongold, n_neg, replace = TRUE)
  )
  ng[, c("gene_a", "gene_b") := list(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
  ng <- unique(ng[gene_a != gene_b])
  ng[, label := 0L]

  pairs <- rbind(gg, ng)

  # Score each pair per layer (binary: in-layer or not)
  for (ln in roc_layers_single) {
    sc <- layer_dt[[ln]][pairs, on = c("gene_a", "gene_b")]
    sc[is.na(in_layer), in_layer := 0L]
    # Skip if degenerate (all-zero or all-one)
    if (length(unique(sc$in_layer)) < 2) next
    r <- suppressMessages(pROC::roc(sc$label, sc$in_layer,
                                    levels = c(0, 1), direction = "<"))
    roc_rows[[length(roc_rows) + 1L]] <- data.table(
      gold_set = gn, scheme = ln, auroc = as.numeric(auc(r)),
      n_pos = sum(sc$label == 1), n_neg = sum(sc$label == 0)
    )
    # PR curve
    ord <- order(-sc$in_layer, runif(nrow(sc)))
    lab_s <- sc$label[ord]
    tp <- cumsum(lab_s == 1); fp <- cumsum(lab_s == 0)
    prec <- tp / (tp + fp); rec <- tp / sum(lab_s == 1)
    # AUPRC via trapezoidal
    auprc <- sum(diff(c(0, rec)) * prec, na.rm = TRUE)
    pr_rows[[length(pr_rows) + 1L]] <- data.table(
      gold_set = gn, scheme = ln, auprc = auprc
    )
  }

  # Composite scoring
  if (has_composite) {
    sc <- comp_dt[pairs, on = c("gene_a", "gene_b")]
    sc[is.na(p_composite), p_composite := 0]
    if (length(unique(sc$p_composite)) >= 2) {
      r <- suppressMessages(pROC::roc(sc$label, sc$p_composite,
                                      levels = c(0, 1), direction = "<"))
      roc_rows[[length(roc_rows) + 1L]] <- data.table(
        gold_set = gn, scheme = "composite_p", auroc = as.numeric(auc(r)),
        n_pos = sum(sc$label == 1), n_neg = sum(sc$label == 0)
      )
      ord <- order(-sc$p_composite, runif(nrow(sc)))
      lab_s <- sc$label[ord]
      tp <- cumsum(lab_s == 1); fp <- cumsum(lab_s == 0)
      prec <- tp / (tp + fp); rec <- tp / sum(lab_s == 1)
      auprc <- sum(diff(c(0, rec)) * prec, na.rm = TRUE)
      pr_rows[[length(pr_rows) + 1L]] <- data.table(
        gold_set = gn, scheme = "composite_p", auprc = auprc
      )
    }
  }
}

roc_tbl <- rbindlist(roc_rows)
pr_tbl  <- rbindlist(pr_rows)
fwrite(roc_tbl, file.path(DATDIR, "bench_roc.csv"))
fwrite(pr_tbl,  file.path(DATDIR, "bench_pr.csv"))

# Plot: AUROC composite vs best single layer per gold set
if (nrow(roc_tbl) > 0) {
  roc_tbl[, is_composite := scheme == "composite_p"]
  roc_tbl[, scheme := factor(scheme,
                             levels = c(roc_layers_single, "composite_p"))]

  p_roc <- ggplot(roc_tbl, aes(x = scheme, y = auroc, fill = is_composite)) +
    geom_col() +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.3, size = 2.8) +
    facet_wrap(~ gold_set, scales = "free_x") +
    scale_fill_manual(values = c(`FALSE` = "#8c96c6", `TRUE` = "#b2182b"),
                      guide = "none") +
    labs(
      title = "Gold-pair recovery: AUROC by scoring scheme",
      subtitle = "Composite p (red) vs per-layer edge presence (grey) on gold-gold vs gold-nongold pairs",
      x = "Scheme", y = "AUROC"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(OUTDIR, "bench_ROC_composite_vs_layers.pdf"),
         p_roc, width = 10, height = 6)
  cat("Wrote: ", file.path(OUTDIR, "bench_ROC_composite_vs_layers.pdf"), "\n",
      sep = "")
}

# ---------------------------------------------------------------------------
# Layer ranking summary
# ---------------------------------------------------------------------------
rank_tbl <- bench[!is.na(fold_over_null) & is.finite(fold_over_null),
                  .(mean_fold = mean(fold_over_null),
                    median_fold = median(fold_over_null),
                    mean_p = mean(p_perm, na.rm = TRUE)),
                  by = layer][order(-mean_fold)]
fwrite(rank_tbl, file.path(DATDIR, "layer_ranking_v2.csv"))

cat("\n=== Layer ranking by mean fold-over-null (degree-preserving) ===\n")
print(rank_tbl)

cat("\nTotal elapsed: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    " mins\n", sep = "")
cat("All outputs -> ", OUTDIR, "\n", sep = "")
