#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --job-name=net_252_coexpr
#SBATCH --output=logs/net_252_coexpr_%j.out
#SBATCH --error=logs/net_252_coexpr_%j.err
# ===========================================================================
# Script 252: Build co-expression edges from merged logCPM matrix
# ===========================================================================
# Purpose: Layer 2 of the 10-layer Bayesian multiplex gene network.
#          Compute pairwise Pearson correlations on the merged logCPM
#          expression matrix (10-cohort integration), keep the top 0.1%
#          by |r|, and output edges with p-values.
#
# Input:
#   - RNA-seq/results/network/network_nodes.csv  (node set V from Script 250)
#   - RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds
#
# Output:
#   - RNA-seq/results/network/edges_coexpr.csv
#     Columns: gene_a, gene_b, raw_score (|Pearson r|), metadata (p-value, n_samples)
#
# Environment: rnaseq (R 4.4+, data.table, arrow, edgeR)
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  # library(arrow)  # using fwrite instead
  library(edgeR)
})

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NODE_PATH <- file.path(BASE, "RNA-seq/results/network/network_nodes.csv")
DGE_PATH  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
OUTDIR    <- file.path(BASE, "RNA-seq/results/network")
OUT_CSV <- file.path(OUTDIR, "edges_coexpr.csv")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 252: Co-expression Edges ===\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---------------------------------------------------------------------------
# 1. Load node set V
# ---------------------------------------------------------------------------
cat("--- Loading node set ---\n")
nodes <- fread(NODE_PATH)
cat("Node set V:", nrow(nodes), "genes\n")

# Strip version suffixes from ensembl_id for matching
nodes[, gene := sub("\\.\\d+$", "", ensembl_id)]
node_genes <- unique(nodes$gene)
cat("Unique unversioned Ensembl IDs:", length(node_genes), "\n\n")

# ---------------------------------------------------------------------------
# 2. Load expression matrix (merged DGEList -> logCPM)
# ---------------------------------------------------------------------------
cat("--- Loading expression data ---\n")
t_load <- Sys.time()
dge <- readRDS(DGE_PATH)
cat("DGEList:", nrow(dge), "genes x", ncol(dge), "samples\n")

logcpm <- cpm(dge, log = TRUE, prior.count = 1)
n_samples <- ncol(logcpm)
cat("logCPM matrix:", nrow(logcpm), "x", n_samples, "\n")
cat("Expression load time:", round(difftime(Sys.time(), t_load, units = "secs"), 1), "s\n\n")

# Free DGEList — only need logcpm from here
rm(dge); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 3. Subset to genes in V
# ---------------------------------------------------------------------------
cat("--- Subsetting to node set ---\n")

# Expression matrix rownames are versioned Ensembl IDs
expr_genes_unversioned <- sub("\\.\\d+$", "", rownames(logcpm))
names(expr_genes_unversioned) <- rownames(logcpm)

# Find versioned rownames that map to our node genes
keep_mask <- expr_genes_unversioned %in% node_genes
logcpm_sub <- logcpm[keep_mask, ]

# Rename rows to unversioned IDs (handle possible duplicates by keeping first)
rownames(logcpm_sub) <- expr_genes_unversioned[rownames(logcpm_sub)]
dup_mask <- duplicated(rownames(logcpm_sub))
if (any(dup_mask)) {
  cat("Dropping", sum(dup_mask), "duplicate unversioned gene IDs\n")
  logcpm_sub <- logcpm_sub[!dup_mask, ]
}

cat("Genes in V AND expression matrix:", nrow(logcpm_sub), "of", length(node_genes), "node genes\n")
cat("Samples:", ncol(logcpm_sub), "\n\n")

# Free full matrix
rm(logcpm); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 4. Compute pairwise Pearson correlation
# ---------------------------------------------------------------------------
cat("--- Computing correlation matrix ---\n")
t_cor <- Sys.time()

n_genes <- nrow(logcpm_sub)
cat("Matrix dimensions:", n_genes, "x", n_genes, "\n")
cat("Estimated memory: ~", round(n_genes^2 * 8 / 1e9, 2), "GB\n")

cor_mat <- cor(t(logcpm_sub), method = "pearson")

cat("Correlation computed in:", round(difftime(Sys.time(), t_cor, units = "secs"), 1), "s\n\n")

# ---------------------------------------------------------------------------
# 5. Determine top 0.1% threshold on |r|
# ---------------------------------------------------------------------------
cat("--- Finding top 0.1% threshold ---\n")

# Zero out diagonal and lower triangle to get unique pairs only
cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA

# Compute threshold from upper triangle values
upper_vals <- cor_mat[upper.tri(cor_mat, diag = FALSE)]
upper_vals <- upper_vals[!is.na(upper_vals)]
n_pairs_total <- length(upper_vals)
cat("Total unique gene pairs:", format(n_pairs_total, big.mark = ","), "\n")

thresh <- quantile(abs(upper_vals), probs = 0.999, na.rm = TRUE)
cat("|r| threshold at top 0.1%:", round(thresh, 4), "\n")

# Free the flat vector
rm(upper_vals); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 6. Extract edges above threshold
# ---------------------------------------------------------------------------
cat("\n--- Extracting edges ---\n")
t_extract <- Sys.time()

idx <- which(abs(cor_mat) >= thresh, arr.ind = TRUE)
cat("Edges above threshold:", format(nrow(idx), big.mark = ","), "\n")

gene_names <- rownames(logcpm_sub)

edges <- data.table(
  gene_a    = gene_names[idx[, 1]],
  gene_b    = gene_names[idx[, 2]],
  raw_score = abs(cor_mat[idx]),
  r_signed  = cor_mat[idx]
)

# Free correlation matrix
rm(cor_mat, idx); gc(verbose = FALSE)

# Ensure gene_a < gene_b alphabetically (canonical ordering), no self-loops
edges <- edges[gene_a != gene_b]
swap <- edges$gene_a > edges$gene_b
if (any(swap)) {
  tmp_a <- edges$gene_a[swap]
  edges[swap, gene_a := gene_b]
  edges[swap, gene_b := tmp_a]
  rm(tmp_a)
}

# Remove any duplicate pairs that might arise from swapping
edges <- unique(edges, by = c("gene_a", "gene_b"))

cat("Edges after dedup + no self-loops:", format(nrow(edges), big.mark = ","), "\n")
cat("Edge extraction time:", round(difftime(Sys.time(), t_extract, units = "secs"), 1), "s\n\n")

# ---------------------------------------------------------------------------
# 7. Compute p-values via t-distribution formula
# ---------------------------------------------------------------------------
cat("--- Computing p-values ---\n")
# t = r * sqrt((n-2) / (1-r^2)), df = n-2
# p = 2 * pt(-abs(t), df = n-2)

n <- n_samples
df <- n - 2L
edges[, t_stat := r_signed * sqrt(df / (1 - r_signed^2))]
edges[, pvalue := 2 * pt(-abs(t_stat), df = df)]

cat("n_samples:", n, "\n")
cat("p-value range: [", format(min(edges$pvalue), digits = 3),
    ",", format(max(edges$pvalue), digits = 3), "]\n\n")

# ---------------------------------------------------------------------------
# 8. Format output and save as parquet
# ---------------------------------------------------------------------------
cat("--- Saving parquet ---\n")

# Construct metadata string: "p={pvalue};n={n_samples}"
edges[, metadata := paste0("p=", signif(pvalue, 4), ";n=", n)]

out <- edges[, .(gene_a, gene_b, raw_score, metadata)]
setorder(out, -raw_score)

fwrite(out, OUT_CSV)
cat("Saved:", OUT_CSV, "\n")
cat("File size:", round(file.size(OUT_CSV) / 1e6, 1), "MB\n\n")

# ---------------------------------------------------------------------------
# 9. Summary statistics
# ---------------------------------------------------------------------------
cat("=== Summary ===\n")
cat("Total edges:", format(nrow(out), big.mark = ","), "\n")
cat("|r| threshold (top 0.1%):", round(thresh, 4), "\n")
cat("Mean |r|:", round(mean(out$raw_score), 4), "\n")
cat("Median |r|:", round(median(out$raw_score), 4), "\n")
cat("Max |r|:", round(max(out$raw_score), 4), "\n")
cat("n_samples:", n, "\n")

# Top 10 highest-degree genes
cat("\n--- Top 10 highest-degree genes ---\n")
degree_a <- edges[, .N, by = gene_a]
degree_b <- edges[, .N, by = gene_b]
setnames(degree_a, c("gene", "deg_a"))
setnames(degree_b, c("gene", "deg_b"))
degree <- merge(degree_a, degree_b, by = "gene", all = TRUE)
degree[is.na(deg_a), deg_a := 0L]
degree[is.na(deg_b), deg_b := 0L]
degree[, degree := deg_a + deg_b]

# Add symbol from node table
sym_map <- nodes[, .(gene, human_symbol)][!duplicated(gene)]
degree <- merge(degree, sym_map, by = "gene", all.x = TRUE)

setorder(degree, -degree)
cat(sprintf("%-18s %-12s %s\n", "ensembl_id", "symbol", "degree"))
cat(strrep("-", 42), "\n")
for (i in seq_len(min(10, nrow(degree)))) {
  cat(sprintf("%-18s %-12s %d\n",
              degree$gene[i],
              ifelse(is.na(degree$human_symbol[i]), ".", degree$human_symbol[i]),
              degree$degree[i]))
}

t1 <- Sys.time()
cat("\nElapsed:", round(difftime(t1, t0, units = "secs"), 1), "s\n")
cat("Done.\n")
