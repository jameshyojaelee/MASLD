#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --job-name=net_291_benchmark
#SBATCH --output=logs/net_291_benchmark_%j.out
#SBATCH --error=logs/net_291_benchmark_%j.err
#
# 291_edge_quality_report.R
# MASLD Network v2 (Architecture C) -- Phase 3: per-edge-type fold-over-null
# benchmark against MASLD gold standards.
#
# Architecture C has NO composite score. Instead, we benchmark each edge type
# independently against gold standards using Maslov-Sneppen degree-preserving
# nulls (igraph::rewire + keeping_degseq). Results are reported as a
# supplementary transparency table -- no hidden integration.
#
# Pre-registration Criterion 3 (docs/network_c_prereg.md):
#   d_coloc fold-over-null on `druggable_druggable` gold pairs must be >= 3x
#   (p_perm < 0.01).
#
# Environment: micromamba activate rnaseq
#
# Inputs:
#   RNA-seq/results/network/edges_ppi.csv          (STRING; filter combined_score >= 0.7)
#   RNA-seq/results/network/edges_d_f2_loco.csv    (or edges_d_f2.csv fallback)
#   RNA-seq/results/network/edges_d_coloc.csv
#   RNA-seq/results/network/edges_d_lr.csv
#   RNA-seq/results/network/edges_d_cerna.csv
#   RNA-seq/results/network/network_nodes.csv
#   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv   (for dgidb_druggable)
#
# Outputs:
#   figures/misc/network_v2_benchmark/data/edge_quality_per_type.csv
#   figures/misc/network_v2_benchmark/edge_quality_heatmap.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(parallel)
})

t0 <- Sys.time()

# ---------------------------------------------------------------------------
# Paths & config
# ---------------------------------------------------------------------------
BASE   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
NETDIR <- file.path(BASE, "RNA-seq/results/network")
OUTDIR <- file.path(BASE, "figures/misc/network_v2_benchmark")
DATDIR <- file.path(OUTDIR, "data")
dir.create(DATDIR, recursive = TRUE, showWarnings = FALSE)

N_REWIRE  <- 1000
SWAP_MULT <- 10
N_CORES   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
if (is.na(N_CORES) || N_CORES < 1) N_CORES <- 8
MC_CORES  <- min(8L, N_CORES)
SEED      <- 42
set.seed(SEED)

PPI_MIN_SCORE <- 700L   # STRING combined_score >= 0.7 (integer scale 0-1000)

cat("======================================================================\n")
cat("291 EDGE QUALITY REPORT (Architecture C, per-edge-type benchmark)\n")
cat("  Cores available: ", N_CORES, "  | mclapply cores: ", MC_CORES, "\n", sep = "")
cat("  N_REWIRE: ", N_REWIRE, "  | SWAP_MULT: ", SWAP_MULT, "\n", sep = "")
cat("======================================================================\n\n")

# ---------------------------------------------------------------------------
# Load nodes
# ---------------------------------------------------------------------------
nodes <- fread(file.path(NETDIR, "network_nodes.csv"))
all_genes <- unique(nodes$human_symbol)
cat("Nodes in V: ", length(all_genes), "\n\n", sep = "")

# ---------------------------------------------------------------------------
# Load edge files. Architecture C edge types:
#   string_ge700 (backbone), d_f2, d_coloc, d_lr, d_cerna
# ---------------------------------------------------------------------------
read_edges <- function(path, score_col = NULL, score_min = NULL) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  if (!is.null(score_col) && score_col %in% names(dt) && !is.null(score_min)) {
    dt <- dt[get(score_col) >= score_min]
  }
  if (!all(c("gene_a", "gene_b") %in% names(dt))) {
    stop("edge file missing gene_a/gene_b: ", path)
  }
  dt <- dt[, .(gene_a, gene_b)]
  dt[, c("gene_a", "gene_b") := list(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
  dt <- unique(dt[gene_a != gene_b & gene_a %in% all_genes & gene_b %in% all_genes])
  dt
}

edge_files <- list()

# STRING PPI backbone (filter to combined_score >= 700)
ppi_path <- file.path(NETDIR, "edges_ppi.csv")
if (file.exists(ppi_path)) {
  ppi_raw <- fread(ppi_path)
  score_col <- intersect(c("combined_score", "score", "string_score"), names(ppi_raw))
  if (length(score_col) > 0) {
    ppi_raw <- ppi_raw[get(score_col[1]) >= PPI_MIN_SCORE]
  }
  ppi_raw <- ppi_raw[, .(gene_a, gene_b)]
  ppi_raw[, c("gene_a", "gene_b") := list(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
  ppi_raw <- unique(ppi_raw[gene_a != gene_b &
                            gene_a %in% all_genes & gene_b %in% all_genes])
  edge_files[["string_ge700"]] <- ppi_raw
}

# d_f2 (prefer LOCO if present)
f2_loco <- file.path(NETDIR, "edges_d_f2_loco.csv")
f2_base <- file.path(NETDIR, "edges_d_f2.csv")
f2_use  <- if (file.exists(f2_loco)) f2_loco else f2_base
if (file.exists(f2_use)) {
  edge_files[["d_f2"]] <- read_edges(f2_use)
  cat("  d_f2 source: ", basename(f2_use), "\n", sep = "")
}

for (ln in c("d_coloc", "d_lr", "d_cerna")) {
  e <- read_edges(file.path(NETDIR, paste0("edges_", ln, ".csv")))
  if (!is.null(e)) edge_files[[ln]] <- e
}

cat("\nEdge types loaded:\n")
for (ln in names(edge_files)) {
  cat(sprintf("  %-14s : %7d edges\n", ln, nrow(edge_files[[ln]])))
}
cat("\n")

# ---------------------------------------------------------------------------
# Gold standards
# ---------------------------------------------------------------------------
pc_path <- file.path(BASE, "results/library/positive_control.csv")
pc <- fread(pc_path)
clin_col <- intersect(c("Clinical_Drug", "clinical_drug"), names(pc))
if (length(clin_col) > 0) {
  sym_col <- intersect(c("Gene symbol", "gene_symbol", "symbol", "gene"), names(pc))
  sym_col <- if (length(sym_col) > 0) sym_col[1] else names(pc)[1]
  gold_clinical <- pc[!is.na(get(clin_col[1])) & get(clin_col[1]) != "", get(sym_col)]
} else {
  gold_clinical <- character(0)
}
gold_clinical <- intersect(gold_clinical, all_genes)

govaere25_published <- c(
  "AKR1B10", "SOX4", "TTC39A", "VIL1", "CLIC6", "ABCB4", "LTBP2", "STMN2",
  "AGMO", "ANGPTL4", "CA3", "CCBE1", "CFB", "CSAD", "FABP4", "FADS1",
  "IGFBP7", "MME", "PANX2", "PLS1", "RP11-1399P15.1", "SDC1", "SRPX",
  "TMEM9", "TPBG"
)
gold_govaere <- intersect(govaere25_published, all_genes)

resm_pathway <- c("THRB", "PCSK9", "PPARA", "APOA1", "CYP7A1")
gold_resm <- intersect(resm_pathway, all_genes)

# druggable_druggable: single gene set; pairs are all (druggable x druggable)
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
atlas <- fread(atlas_path)
sym_col_atlas <- intersect(c("gene_symbol", "human_symbol", "symbol", "gene"),
                           names(atlas))
if (length(sym_col_atlas) == 0) stop("atlas missing gene symbol column")
drug_col <- intersect(c("dgidb_druggable", "DGIdb_druggable", "is_druggable"),
                      names(atlas))
if (length(drug_col) == 0) stop("atlas missing dgidb_druggable column")
gold_drug <- atlas[as.logical(get(drug_col[1])) %in% TRUE |
                   get(drug_col[1]) %in% c(1, "1", "TRUE", "True", "yes", "Yes"),
                   get(sym_col_atlas[1])]
gold_drug <- intersect(unique(gold_drug), all_genes)

gold_sets <- list(
  clinical_drug       = gold_clinical,
  govaere25           = gold_govaere,
  resmetirom_pathway  = gold_resm,
  druggable_druggable = gold_drug
)

cat("Gold-set sizes (in V):\n")
for (gn in names(gold_sets)) {
  cat(sprintf("  %-22s : %5d\n", gn, length(gold_sets[[gn]])))
}
cat("\n")

# ---------------------------------------------------------------------------
# Helpers (Maslov-Sneppen rewiring; pattern from 272b)
# ---------------------------------------------------------------------------
gs_connectivity <- function(el, gold_genes) {
  if (length(gold_genes) < 2 || nrow(el) == 0) return(NA_real_)
  n_pairs_total <- choose(length(gold_genes), 2)
  gs <- el[gene_a %in% gold_genes & gene_b %in% gold_genes]
  nrow(gs) / n_pairs_total
}

build_graph <- function(el) {
  if (nrow(el) == 0) return(NULL)
  g <- graph_from_data_frame(el[, .(gene_a, gene_b)], directed = FALSE,
                             vertices = data.frame(name = all_genes))
  simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

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

# Per-edge-type worker: computes null fracs for all gold sets given a graph
benchmark_edge_type <- function(ln, el, gold_sets, seed_base) {
  g <- build_graph(el)
  n_edges <- ecount(g)
  n_swaps <- SWAP_MULT * n_edges
  rows <- list()
  for (gn in names(gold_sets)) {
    gold <- gold_sets[[gn]]
    if (length(gold) < 2) next
    observed <- gs_connectivity(el, gold)
    t_null <- Sys.time()
    null_fracs <- vapply(seq_len(N_REWIRE), function(i) {
      rewire_one(g, n_swaps, gold, seed_i = seed_base + i)
    }, numeric(1))
    null_time <- as.numeric(difftime(Sys.time(), t_null, units = "secs"))

    mean_null <- mean(null_fracs, na.rm = TRUE)
    sd_null   <- sd(null_fracs, na.rm = TRUE)
    fold <- if (!is.na(mean_null) && mean_null > 0) observed / mean_null else NA_real_
    n_ge <- sum(null_fracs >= observed, na.rm = TRUE)
    p_perm <- (n_ge + 1) / (sum(!is.na(null_fracs)) + 1)

    fold_dist <- observed / null_fracs
    fold_dist <- fold_dist[is.finite(fold_dist)]
    ci_lo <- if (length(fold_dist) > 10) quantile(fold_dist, 0.025, na.rm = TRUE) else NA_real_
    ci_hi <- if (length(fold_dist) > 10) quantile(fold_dist, 0.975, na.rm = TRUE) else NA_real_

    rows[[length(rows) + 1L]] <- data.table(
      edge_type      = ln,
      gold_set       = gn,
      n_edges        = n_edges,
      n_gold         = length(gold),
      n_gold_pairs   = choose(length(gold), 2),
      observed_frac  = observed,
      null_mean      = mean_null,
      null_sd        = sd_null,
      null_95ci_lo   = ci_lo,
      null_95ci_hi   = ci_hi,
      fold_over_null = fold,
      p_perm         = p_perm,
      n_null         = sum(!is.na(null_fracs)),
      null_seconds   = null_time
    )
  }
  rbindlist(rows)
}

# ---------------------------------------------------------------------------
# Parallelize across edge types (each type is independent)
# ---------------------------------------------------------------------------
edge_type_names <- names(edge_files)
cat("Running benchmark across ", length(edge_type_names), " edge types (mc.cores=",
    MC_CORES, ")\n\n", sep = "")

results_list <- mclapply(seq_along(edge_type_names), function(i) {
  ln <- edge_type_names[i]
  el <- edge_files[[ln]]
  cat(">> starting: ", ln, " (", nrow(el), " edges)\n", sep = "")
  t_layer <- Sys.time()
  out <- benchmark_edge_type(ln, el, gold_sets,
                             seed_base = SEED + 1000L * i)
  cat("<< done: ", ln, " in ",
      round(as.numeric(difftime(Sys.time(), t_layer, units = "mins")), 2),
      " min\n", sep = "")
  out
}, mc.cores = MC_CORES, mc.preschedule = FALSE)

bench <- rbindlist(results_list)

out_csv <- file.path(DATDIR, "edge_quality_per_type.csv")
fwrite(bench, out_csv)
cat("\nWrote: ", out_csv, "\n\n", sep = "")

# Console summary
cat("=== Per-edge-type fold-over-null ===\n")
print(bench[, .(edge_type, gold_set, n_edges, n_gold_pairs,
                observed_frac, null_mean, fold_over_null,
                null_95ci_lo, null_95ci_hi, p_perm)])
cat("\n")

# ---------------------------------------------------------------------------
# Figure: heatmap of log10(fold_over_null)
# ---------------------------------------------------------------------------
bench_plot <- copy(bench)
bench_plot <- bench_plot[!is.na(fold_over_null) & is.finite(fold_over_null) &
                           fold_over_null > 0]

edge_order <- intersect(c("string_ge700", "d_f2", "d_coloc", "d_lr", "d_cerna"),
                        unique(bench_plot$edge_type))
gs_order <- intersect(c("clinical_drug", "govaere25", "resmetirom_pathway",
                        "druggable_druggable"),
                      unique(bench_plot$gold_set))
bench_plot[, edge_type := factor(edge_type, levels = edge_order)]
bench_plot[, gold_set  := factor(gold_set,  levels = gs_order)]
bench_plot[, log10_fold := log10(pmax(fold_over_null, 1e-3))]
bench_plot[, label := sprintf("%.1fx\n[%.1f-%.1f]",
                              fold_over_null, null_95ci_lo, null_95ci_hi)]
bench_plot[!is.finite(null_95ci_lo) | !is.finite(null_95ci_hi),
           label := sprintf("%.1fx", fold_over_null)]

p_heat <- ggplot(bench_plot, aes(x = gold_set, y = edge_type, fill = log10_fold)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 2.8, color = "black", lineheight = 0.9) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, name = "log10 fold\n(deg-pres null)") +
  labs(
    title = "Architecture C -- per-edge-type fold-over-null",
    subtitle = sprintf("N=%d Maslov-Sneppen rewirings per edge type; 95%% CI in brackets",
                       N_REWIRE),
    x = "Gold-standard gene set", y = "Edge type"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid  = element_blank())

pdf_path <- file.path(OUTDIR, "edge_quality_heatmap.pdf")
ggsave(pdf_path, p_heat, width = 9, height = 5)
cat("Wrote: ", pdf_path, "\n\n", sep = "")

# ---------------------------------------------------------------------------
# Criterion 3 check: d_coloc on druggable_druggable fold >= 3x AND p_perm < 0.01
# ---------------------------------------------------------------------------
c3 <- bench[edge_type == "d_coloc" & gold_set == "druggable_druggable"]
cat("======================================================================\n")
if (nrow(c3) == 0) {
  cat("CRITERION 3: UNKNOWN -- d_coloc x druggable_druggable row missing\n")
} else {
  fold_c3 <- c3$fold_over_null[1]
  p_c3    <- c3$p_perm[1]
  pass    <- !is.na(fold_c3) && !is.na(p_c3) && fold_c3 >= 3 && p_c3 < 0.01
  status  <- if (pass) "PASS" else "FAIL"
  cat(sprintf("CRITERION 3: %s  |  d_coloc x druggable_druggable fold=%.2fx  p_perm=%.2e\n",
              status, fold_c3, p_c3))
  cat(sprintf("  Threshold: fold >= 3x AND p_perm < 0.01\n"))
}
cat("======================================================================\n\n")

cat("Total elapsed: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    " mins\n", sep = "")
cat("All outputs -> ", OUTDIR, "\n", sep = "")
