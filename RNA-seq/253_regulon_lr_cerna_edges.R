#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_253_reg_lr_cerna
#SBATCH --output=logs/net_253_reg_lr_cerna_%j.out
#SBATCH --error=logs/net_253_reg_lr_cerna_%j.err
# 253_regulon_lr_cerna_edges.R — Build three edge types for the Bayesian multiplex network:
#   Layer 3: TF regulon co-membership (Jaccard of shared regulators)
#   Layer 4: Ligand-receptor interactions (LIANA differential scores)
#   Layer 9: ceRNA shared miRNA edges (lncRNA-mRNA pairs)

suppressPackageStartupMessages({
  library(data.table)
  # library(arrow)  # using fwrite instead
})

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR     <- file.path(BASE, "RNA-seq/results/network")
NODE_PATH  <- file.path(OUTDIR, "network_nodes.csv")
REG_PATH   <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")
LIANA_PATH <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv")
CERNA_PATH <- file.path(BASE, "RNA-seq/results/ncrna/cerna_network_full.csv")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 253: Regulon / LR / ceRNA Edge Construction ===\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---------- Load node set V ----------
stopifnot(file.exists(NODE_PATH))
nodes <- fread(NODE_PATH)
V <- unique(nodes$human_symbol)
cat("Node set V:", length(V), "genes\n\n")

# ==========================================================================
#  LAYER 3 — TF Regulon co-membership
# ==========================================================================
cat("--- Layer 3: TF Regulon Co-membership ---\n")

reg <- fread(REG_PATH)
# Columns: tf_name, target_gene, regulon_id, n_target_genes, target_genes,
#           n_enhancers, mean_activity_masld, mean_activity_normal,
#           regulon_activity_diff, activity_pval, activity_padj

cat("  Regulon file:", nrow(reg), "rows,", uniqueN(reg$tf_name), "TFs\n")

# Build gene -> set-of-regulators mapping (including TF self-membership)
gene_to_tfs <- reg[, .(tfs = list(unique(tf_name))), by = target_gene]
setnames(gene_to_tfs, "target_gene", "gene")

# Restrict to genes in V
gene_to_tfs <- gene_to_tfs[gene %in% V]
cat("  Genes in V with regulon membership:", nrow(gene_to_tfs), "\n")

# --- 3a: TF-target edges (TF -> target in same regulon) ---
tf_target_edges <- reg[target_gene %in% V & tf_name %in% V & tf_name != target_gene,
                       .(gene_a = tf_name, gene_b = target_gene,
                         regulon_id, activity_padj)]
# Use 1 - padj as confidence (lower padj = higher confidence); cap at 1.0
tf_target_edges[, raw_score := fifelse(is.na(activity_padj), 1.0,
                                       pmin(1.0, 1.0 - activity_padj))]
# Deduplicate (a TF may appear as its own target via self-regulation row)
tf_target_edges <- unique(tf_target_edges, by = c("gene_a", "gene_b"))
tf_target_edges[, metadata := paste0("tf_target|regulon=", regulon_id)]

# --- 3b: Co-regulated gene pairs (Jaccard of shared TF regulons) ---
if (nrow(gene_to_tfs) >= 2) {
  genes_with_reg <- gene_to_tfs$gene
  n_genes <- length(genes_with_reg)

  # Pre-compute TF sets as character vectors for fast Jaccard
  tf_sets <- setNames(gene_to_tfs$tfs, gene_to_tfs$gene)

  # Build all pairs (upper triangle)
  pairs_idx <- combn(seq_len(n_genes), 2)
  n_pairs <- ncol(pairs_idx)

  gene_a_vec <- character(n_pairs)
  gene_b_vec <- character(n_pairs)
  jaccard_vec <- numeric(n_pairs)
  shared_str <- character(n_pairs)
  n_shared_vec <- integer(n_pairs)

  for (k in seq_len(n_pairs)) {
    i <- pairs_idx[1, k]
    j <- pairs_idx[2, k]
    s_i <- tf_sets[[i]]
    s_j <- tf_sets[[j]]
    shared <- intersect(s_i, s_j)
    n_sh <- length(shared)
    if (n_sh > 0) {
      gene_a_vec[k] <- genes_with_reg[i]
      gene_b_vec[k] <- genes_with_reg[j]
      jaccard_vec[k] <- n_sh / length(union(s_i, s_j))
      shared_str[k] <- paste(sort(shared), collapse = ";")
      n_shared_vec[k] <- n_sh
    }
  }

  coreg_edges <- data.table(gene_a = gene_a_vec, gene_b = gene_b_vec,
                            raw_score = jaccard_vec, shared_regulons = shared_str,
                            n_shared = n_shared_vec)
  coreg_edges <- coreg_edges[raw_score > 0]
  coreg_edges[, metadata := paste0("co_regulated|shared=", shared_regulons, "|n=", n_shared)]
} else {
  coreg_edges <- data.table(gene_a = character(), gene_b = character(),
                            raw_score = numeric(), shared_regulons = character(),
                            n_shared = integer(), metadata = character())
}

# Combine TF-target and co-regulation edges
regulon_out <- rbind(
  tf_target_edges[, .(gene_a, gene_b, raw_score, metadata)],
  coreg_edges[, .(gene_a, gene_b, raw_score, metadata)],
  fill = TRUE
)

# Ensure canonical ordering (gene_a < gene_b alphabetically) then deduplicate
regulon_out[, c("ga", "gb") := .(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
regulon_out[, c("gene_a", "gene_b") := .(ga, gb)]
regulon_out[, c("ga", "gb") := NULL]
# Keep the higher raw_score when duplicated
setorder(regulon_out, gene_a, gene_b, -raw_score)
regulon_out <- regulon_out[!duplicated(regulon_out, by = c("gene_a", "gene_b"))]

fwrite(regulon_out, file.path(OUTDIR, "edges_regulon.csv"))

cat("  TF-target edges (both in V):", nrow(tf_target_edges), "\n")
cat("  Co-regulation edges (Jaccard>0):", nrow(coreg_edges), "\n")
cat("  Combined unique edges:", nrow(regulon_out), "\n")
cat("  Mean raw_score:", round(mean(regulon_out$raw_score), 4), "\n")
if (nrow(regulon_out) > 0) {
  top_reg <- regulon_out[order(-raw_score)][1:min(5, .N)]
  cat("  Top edges:\n")
  for (r in seq_len(nrow(top_reg))) {
    cat("    ", top_reg$gene_a[r], " -- ", top_reg$gene_b[r],
        " (score=", round(top_reg$raw_score[r], 4), ")\n")
  }
}
cat("\n")

# ==========================================================================
#  LAYER 4 — Ligand-Receptor interactions
# ==========================================================================
cat("--- Layer 4: Ligand-Receptor Interactions (LIANA) ---\n")

liana <- fread(LIANA_PATH)
# Columns: source, target, ligand_complex, receptor_complex,
#           score_masld, score_control, score_diff

cat("  LIANA file:", nrow(liana), "rows,",
    uniqueN(liana$ligand_complex), "ligands,",
    uniqueN(liana$receptor_complex), "receptors\n")

# Decompose complexes: split on "_" for multi-subunit complexes (e.g., B2M_FCGRT)
# Each subunit becomes a separate gene entry
expand_complex <- function(x) {
  strsplit(x, "_")
}

liana[, ligand_genes := expand_complex(ligand_complex)]
liana[, receptor_genes := expand_complex(receptor_complex)]

# Explode to individual gene pairs
lr_expanded <- liana[, {
  lg <- ligand_genes[[1]]
  rg <- receptor_genes[[1]]
  CJ(lig_gene = lg, rec_gene = rg)
}, by = .(source, target, ligand_complex, receptor_complex,
          score_masld, score_control, score_diff)]

# Filter to genes in V
lr_in_v <- lr_expanded[lig_gene %in% V & rec_gene %in% V & lig_gene != rec_gene]

cat("  Expanded LR pairs in V:", nrow(lr_in_v), "\n")

# Aggregate across cell-type pairs: for each ligand-receptor gene pair,
# take the maximum absolute score_diff and record cell-type contexts
lr_agg <- lr_in_v[, .(
  raw_score     = max(abs(score_diff)),
  n_celltype_pairs = .N,
  cell_type_pairs  = paste(unique(paste0(source, "->", target)), collapse = ";"),
  direction        = fifelse(sum(score_diff) > 0, "MASLD_up", "MASLD_down")
), by = .(lig_gene, rec_gene)]

setnames(lr_agg, c("lig_gene", "rec_gene"), c("gene_a", "gene_b"))

# Canonical ordering
lr_agg[, c("ga", "gb") := .(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
lr_agg[, c("gene_a", "gene_b") := .(ga, gb)]
lr_agg[, c("ga", "gb") := NULL]
setorder(lr_agg, gene_a, gene_b, -raw_score)
lr_agg <- lr_agg[!duplicated(lr_agg, by = c("gene_a", "gene_b"))]

# Build metadata string
lr_agg[, metadata := paste0("lr|celltype_pairs=", cell_type_pairs,
                            "|direction=", direction,
                            "|n_contexts=", n_celltype_pairs)]

lr_out <- lr_agg[, .(gene_a, gene_b, raw_score, metadata)]
fwrite(lr_out, file.path(OUTDIR, "edges_lr.csv"))

cat("  Unique LR edges:", nrow(lr_out), "\n")
cat("  Mean raw_score:", round(mean(lr_out$raw_score), 4), "\n")
if (nrow(lr_out) > 0) {
  top_lr <- lr_out[order(-raw_score)][1:min(5, .N)]
  cat("  Top edges:\n")
  for (r in seq_len(nrow(top_lr))) {
    cat("    ", top_lr$gene_a[r], " -- ", top_lr$gene_b[r],
        " (score=", round(top_lr$raw_score[r], 4), ")\n")
  }
}
cat("\n")

# ==========================================================================
#  LAYER 9 — ceRNA shared miRNA edges
# ==========================================================================
cat("--- Layer 9: ceRNA Shared miRNA Edges ---\n")

cerna <- fread(CERNA_PATH)
# Columns: lncrna, mrna, shared_mirnas, cerna_score

cat("  ceRNA file:", nrow(cerna), "triplets\n")

# Filter to genes in V (both lncRNA and mRNA must be present)
cerna_in_v <- cerna[lncrna %in% V & mrna %in% V]
cat("  Triplets with both genes in V:", nrow(cerna_in_v), "\n")

if (nrow(cerna_in_v) > 0) {
  # Count max possible shared miRNAs across all triplets for normalization
  max_possible <- max(cerna_in_v$cerna_score)

  cerna_out <- cerna_in_v[, .(
    gene_a    = lncrna,
    gene_b    = mrna,
    raw_score = cerna_score / max_possible,
    metadata  = paste0("cerna|shared_mirnas=", shared_mirnas,
                       "|n_shared=", cerna_score)
  )]

  # Canonical ordering
  cerna_out[, c("ga", "gb") := .(pmin(gene_a, gene_b), pmax(gene_a, gene_b))]
  cerna_out[, c("gene_a", "gene_b") := .(ga, gb)]
  cerna_out[, c("ga", "gb") := NULL]
  cerna_out <- unique(cerna_out, by = c("gene_a", "gene_b"))
} else {
  cerna_out <- data.table(gene_a = character(), gene_b = character(),
                          raw_score = numeric(), metadata = character())
}

fwrite(cerna_out, file.path(OUTDIR, "edges_cerna.csv"))

cat("  ceRNA edges:", nrow(cerna_out), "\n")
if (nrow(cerna_out) > 0) {
  cat("  Mean raw_score:", round(mean(cerna_out$raw_score), 4), "\n")
  cat("  Edges:\n")
  for (r in seq_len(nrow(cerna_out))) {
    cat("    ", cerna_out$gene_a[r], " -- ", cerna_out$gene_b[r],
        " (score=", round(cerna_out$raw_score[r], 4),
        ", ", cerna_out$metadata[r], ")\n")
  }
}
cat("\n")

# ==========================================================================
#  SUMMARY
# ==========================================================================
cat("=== Summary ===\n")
cat("  Layer 3 (Regulon):  ", nrow(regulon_out), "edges\n")
cat("  Layer 4 (LR):       ", nrow(lr_out), "edges\n")
cat("  Layer 9 (ceRNA):    ", nrow(cerna_out), "edges\n")
cat("  Total:              ", nrow(regulon_out) + nrow(lr_out) + nrow(cerna_out), "edges\n\n")

# Gene coverage across layers
all_genes <- unique(c(regulon_out$gene_a, regulon_out$gene_b,
                      lr_out$gene_a, lr_out$gene_b,
                      cerna_out$gene_a, cerna_out$gene_b))
cat("  Genes touched by any layer:", length(all_genes), "/", length(V),
    "(", round(100 * length(all_genes) / length(V), 1), "%)\n")

cat("\nOutputs:\n")
cat("  ", file.path(OUTDIR, "edges_regulon.parquet"), "\n")
cat("  ", file.path(OUTDIR, "edges_lr.parquet"), "\n")
cat("  ", file.path(OUTDIR, "edges_cerna.parquet"), "\n")

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat("\nDone in", round(as.numeric(elapsed), 1), "minutes.\n")
