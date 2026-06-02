#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_255_pathway
#SBATCH --output=logs/net_255_pathway_%j.out
#SBATCH --error=logs/net_255_pathway_%j.err
# ===========================================================================
# Script 255: Build pathway co-membership edges from MSigDB
# ===========================================================================
# Purpose: Construct pathway co-membership edge layer (Layer 6) for the
#          Bayesian multiplex network. Two genes are connected if they share
#          at least one MSigDB pathway; edge weight = Jaccard index of shared
#          pathway sets.
#
# Input:
#   - RNA-seq/results/network/network_nodes.csv  (node set with human_symbol)
#   - MSigDB gene sets via msigdbr (v10+):
#       Hallmark (H), KEGG (C2:CP:KEGG_MEDICUS or KEGG_LEGACY),
#       Reactome (C2:CP:REACTOME)
#
# Output:
#   - RNA-seq/results/network/edges_pathway.csv
#     Columns: gene_a, gene_b, raw_score (Jaccard), metadata (JSON string)
#
# Efficiency:
#   Instead of O(N^2) pairwise comparisons, uses an inverted index
#   (pathway -> member genes). For each pathway, enumerate all member gene
#   pairs -- these are the ONLY pairs with Jaccard > 0. Accumulate
#   co-membership counts, then compute Jaccard = intersect / union.
#   Complexity: O(sum |pathway|^2) which is much smaller than O(N^2).
#
# Environment: rnaseq (R 4.4+, data.table, arrow, msigdbr v10+)
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  # library(arrow)  # using fwrite instead
  library(msigdbr)
})

t0 <- Sys.time()
cat("=== Script 255: Pathway Co-Membership Edges ===\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n")

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

node_file   <- file.path(BASE, "RNA-seq/results/network/network_nodes.csv")
out_csv <- file.path(BASE, "RNA-seq/results/network/edges_pathway.csv")

stopifnot(file.exists(node_file))

# ---------------------------------------------------------------------------
# 1. Load network node set
# ---------------------------------------------------------------------------
cat("\n--- Loading network nodes ---\n")
nodes <- fread(node_file)
cat("Node file:", nrow(nodes), "rows x", ncol(nodes), "columns\n")

V <- unique(nodes$human_symbol)
V <- V[!is.na(V) & V != ""]
cat("Unique gene symbols in V:", length(V), "\n")

# ---------------------------------------------------------------------------
# 2. Load MSigDB gene sets: Hallmark + KEGG + Reactome
# ---------------------------------------------------------------------------
cat("\n--- Loading MSigDB gene sets ---\n")
t_msigdb <- Sys.time()

# Hallmark
hallmark_df <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
cat("  Hallmark: ", uniqueN(hallmark_df$gs_name), " gene sets,",
    format(nrow(hallmark_df), big.mark = ","), "gene-set pairs\n")

# KEGG (KEGG_MEDICUS preferred, fall back to KEGG_LEGACY)
kegg_df <- tryCatch(
  as.data.table(msigdbr(species = "Homo sapiens",
                         collection = "C2", subcollection = "CP:KEGG_MEDICUS")),
  error = function(e) {
    cat("  KEGG_MEDICUS unavailable, using KEGG_LEGACY\n")
    as.data.table(msigdbr(species = "Homo sapiens",
                           collection = "C2", subcollection = "CP:KEGG_LEGACY"))
  }
)
cat("  KEGG:     ", uniqueN(kegg_df$gs_name), " gene sets,",
    format(nrow(kegg_df), big.mark = ","), "gene-set pairs\n")

# Reactome
reactome_df <- as.data.table(msigdbr(species = "Homo sapiens",
                                      collection = "C2", subcollection = "CP:REACTOME"))
cat("  Reactome: ", uniqueN(reactome_df$gs_name), " gene sets,",
    format(nrow(reactome_df), big.mark = ","), "gene-set pairs\n")

# Combine
gs_all <- rbindlist(list(
  hallmark_df[, .(gs_name, gene_symbol)],
  kegg_df[,     .(gs_name, gene_symbol)],
  reactome_df[, .(gs_name, gene_symbol)]
))
gs_all <- unique(gs_all)  # deduplicate in case of overlap

# Filter to genes in V
gs_all <- gs_all[gene_symbol %chin% V]

n_pathways <- uniqueN(gs_all$gs_name)
n_genes_in_gs <- uniqueN(gs_all$gene_symbol)
cat(sprintf("\nAfter filtering to V: %d pathways covering %d genes (%d gene-pathway pairs)\n",
            n_pathways, n_genes_in_gs, nrow(gs_all)))
cat("MSigDB loading time:", format(round(difftime(Sys.time(), t_msigdb, units = "secs"), 1)), "\n")

# ---------------------------------------------------------------------------
# 3. Build gene -> pathway set membership
# ---------------------------------------------------------------------------
cat("\n--- Building gene -> pathway membership ---\n")
gene2pw <- gs_all[, .(pathways = list(gs_name)), by = gene_symbol]
setkey(gene2pw, gene_symbol)

# Also compute pathway sizes (for ranking shared pathways by specificity)
pw_sizes <- gs_all[, .(size = uniqueN(gene_symbol)), by = gs_name]
setkey(pw_sizes, gs_name)

# Number of pathways per gene (for union size computation)
n_pw_per_gene <- gs_all[, .(n_pw = uniqueN(gs_name)), by = gene_symbol]
setkey(n_pw_per_gene, gene_symbol)
cat("Pathways per gene -- min:", min(n_pw_per_gene$n_pw),
    " median:", median(n_pw_per_gene$n_pw),
    " max:", max(n_pw_per_gene$n_pw), "\n")

# ---------------------------------------------------------------------------
# 4. Build inverted index: pathway -> gene list, enumerate co-member pairs
# ---------------------------------------------------------------------------
cat("\n--- Building co-membership counts via inverted index ---\n")
t_inv <- Sys.time()

# Inverted index: for each pathway, list of member genes in V
pw2genes <- gs_all[, .(genes = list(sort(unique(gene_symbol)))), by = gs_name]

# Pre-allocate list for co-occurrence pairs
# For each pathway, generate all choose(n,2) pairs and accumulate counts
# Use data.table rbindlist for efficiency
cat("Enumerating pairs from", nrow(pw2genes), "pathways...\n")

pair_chunks <- vector("list", nrow(pw2genes))
n_pairs_total <- 0L

for (i in seq_len(nrow(pw2genes))) {
  g <- pw2genes$genes[[i]]
  ng <- length(g)
  if (ng < 2L) next

  # Generate all pairs (gene_a < gene_b guaranteed because g is sorted)
  idx <- combn(ng, 2L)
  pair_chunks[[i]] <- data.table(
    gene_a = g[idx[1L, ]],
    gene_b = g[idx[2L, ]]
  )
  n_pairs_total <- n_pairs_total + ncol(idx)

  if (i %% 500L == 0L) {
    cat(sprintf("  Processed %d / %d pathways (%s pair records so far)\n",
                i, nrow(pw2genes), format(n_pairs_total, big.mark = ",")))
  }
}

cat(sprintf("Total pair records before aggregation: %s\n",
            format(n_pairs_total, big.mark = ",")))

# Combine all pair chunks
pairs <- rbindlist(pair_chunks, use.names = TRUE)
rm(pair_chunks)
gc(verbose = FALSE)

cat("Inverted index + pair enumeration time:",
    format(round(difftime(Sys.time(), t_inv, units = "secs"), 1)), "\n")

# ---------------------------------------------------------------------------
# 5. Aggregate: count shared pathways per gene pair
# ---------------------------------------------------------------------------
cat("\n--- Aggregating co-membership counts ---\n")
t_agg <- Sys.time()

co_counts <- pairs[, .(n_shared = .N), by = .(gene_a, gene_b)]
rm(pairs)
gc(verbose = FALSE)

cat("Unique gene pairs with >= 1 shared pathway:", format(nrow(co_counts), big.mark = ","), "\n")
cat("Aggregation time:", format(round(difftime(Sys.time(), t_agg, units = "secs"), 1)), "\n")

# ---------------------------------------------------------------------------
# 6. Compute Jaccard = n_shared / (n_pw_a + n_pw_b - n_shared)
# ---------------------------------------------------------------------------
cat("\n--- Computing Jaccard indices ---\n")
t_jac <- Sys.time()

# Lookup pathway counts for each gene
co_counts[n_pw_per_gene, n_pw_a := i.n_pw, on = .(gene_a = gene_symbol)]
co_counts[n_pw_per_gene, n_pw_b := i.n_pw, on = .(gene_b = gene_symbol)]

co_counts[, jaccard := n_shared / (n_pw_a + n_pw_b - n_shared)]

cat("Jaccard computation time:", format(round(difftime(Sys.time(), t_jac, units = "secs"), 1)), "\n")

# ---------------------------------------------------------------------------
# 7. Identify top 3 shared pathways per pair (prefer smaller/more specific)
# ---------------------------------------------------------------------------
cat("\n--- Identifying top shared pathways per pair ---\n")
t_meta <- Sys.time()

# Build a fast lookup: for each gene, its pathway set
gene_pw_list <- setNames(gene2pw$pathways, gene2pw$gene_symbol)

# For efficiency on large edge sets, process in batches
batch_size <- 100000L
n_edges <- nrow(co_counts)
top_pw_vec <- character(n_edges)

for (start in seq(1L, n_edges, by = batch_size)) {
  end <- min(start + batch_size - 1L, n_edges)
  batch_idx <- start:end

  top_pw_vec[batch_idx] <- vapply(batch_idx, function(j) {
    pw_a <- gene_pw_list[[co_counts$gene_a[j]]]
    pw_b <- gene_pw_list[[co_counts$gene_b[j]]]
    shared <- intersect(pw_a, pw_b)
    if (length(shared) == 0L) return("")
    # Rank by pathway size (ascending = more specific first)
    shared_sizes <- pw_sizes[.(shared), size, nomatch = NA_integer_]
    shared_sorted <- shared[order(shared_sizes, na.last = TRUE)]
    paste(head(shared_sorted, 3L), collapse = "; ")
  }, character(1L))

  if (start %% 500000L == 1L || end == n_edges) {
    cat(sprintf("  Metadata: %s / %s edges\n",
                format(end, big.mark = ","), format(n_edges, big.mark = ",")))
  }
}

co_counts[, top_shared_pathways := top_pw_vec]
rm(top_pw_vec)

cat("Metadata time:", format(round(difftime(Sys.time(), t_meta, units = "mins"), 1)), "\n")

# ---------------------------------------------------------------------------
# 8. Build final output table
# ---------------------------------------------------------------------------
cat("\n--- Building output table ---\n")

edges <- co_counts[, .(
  gene_a    = gene_a,
  gene_b    = gene_b,
  raw_score = jaccard,
  metadata  = sprintf('{"n_shared_pathways":%d,"top_shared_pathway":"%s"}',
                       n_shared, top_shared_pathways)
)]

# ---------------------------------------------------------------------------
# 9. Save as parquet
# ---------------------------------------------------------------------------
cat("\n--- Saving parquet ---\n")
fwrite(edges, out_csv)
cat("Saved:", out_csv, "\n")
cat("File size:", sprintf("%.1f MB", file.size(out_csv) / 1e6), "\n")

# ---------------------------------------------------------------------------
# 10. Summary statistics
# ---------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Pathways used:", n_pathways, "\n")
cat("  Hallmark: ", uniqueN(hallmark_df$gs_name), "\n")
cat("  KEGG:     ", uniqueN(kegg_df$gs_name), "\n")
cat("  Reactome: ", uniqueN(reactome_df$gs_name), "\n")
cat("Genes in V with pathway annotation:", n_genes_in_gs,
    sprintf("(%.1f%% of %d)\n", 100 * n_genes_in_gs / length(V), length(V)))
cat("Total edges (Jaccard > 0):", format(nrow(edges), big.mark = ","), "\n")

cat("\nJaccard distribution:\n")
cat(sprintf("  Min:    %.4f\n", min(edges$raw_score)))
cat(sprintf("  Q25:    %.4f\n", quantile(edges$raw_score, 0.25)))
cat(sprintf("  Median: %.4f\n", median(edges$raw_score)))
cat(sprintf("  Mean:   %.4f\n", mean(edges$raw_score)))
cat(sprintf("  Q75:    %.4f\n", quantile(edges$raw_score, 0.75)))
cat(sprintf("  Max:    %.4f\n", max(edges$raw_score)))

# Degree distribution
degree <- rbind(
  edges[, .(gene = gene_a)],
  edges[, .(gene = gene_b)]
)[, .N, by = gene][order(-N)]

cat("\nTop 10 highest-degree genes (pathway co-membership):\n")
print(degree[1:min(10, nrow(degree))])

cat(sprintf("\nDegree distribution: min=%d  Q25=%d  median=%d  Q75=%d  max=%d\n",
            min(degree$N), quantile(degree$N, 0.25),
            median(degree$N), quantile(degree$N, 0.75),
            max(degree$N)))

cat("Unique genes in edges:", nrow(degree),
    sprintf("(%.1f%% of V)\n", 100 * nrow(degree) / length(V)))

# Shared pathway count distribution
cat("\nShared pathway count distribution:\n")
cat(sprintf("  Min:    %d\n", min(co_counts$n_shared)))
cat(sprintf("  Median: %d\n", median(co_counts$n_shared)))
cat(sprintf("  Mean:   %.1f\n", mean(co_counts$n_shared)))
cat(sprintf("  Max:    %d\n", max(co_counts$n_shared)))

t1 <- Sys.time()
cat(sprintf("\nTotal elapsed: %s\n", format(round(difftime(t1, t0, units = "mins"), 1))))
cat("Done:", format(t1, "%Y-%m-%d %H:%M:%S"), "\n")
