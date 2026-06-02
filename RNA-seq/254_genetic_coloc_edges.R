#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_254_genetic
#SBATCH --output=logs/net_254_genetic_%j.out
#SBATCH --error=logs/net_254_genetic_%j.err
# ===========================================================================
# Script 254: Build genetic co-localization edges from SuSiE-COLOC results
# ===========================================================================
# Purpose: Layer 5 of the Bayesian multiplex gene network.
#          Two genes are connected if they both colocalize (PP.H4 > 0.5)
#          at the SAME GWAS locus. Edge weight = min(PP.H4_gene_a, PP.H4_gene_b)
#          -- the bottleneck confidence.
#
# Locus definition: genes share a locus if they are tested in the same GWAS,
#   on the same chromosome, and their top_snp positions are within 1 Mb.
#   Single-linkage clustering groups nearby top_snps into contiguous loci.
#
# Input:
#   - RNA-seq/results/network/network_nodes.csv              (node set V)
#   - GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv
#     (322K gene-GWAS pairs, columns: gene, ensembl, chr, gwas_name,
#      PP.H4.abf, top_snp)
#
# Output:
#   - RNA-seq/results/network/edges_genetic.csv
#     Columns: gene_a, gene_b, raw_score (min PP.H4), metadata (shared_gwas, locus_id)
#
# Environment: rnaseq (R 4.4+, data.table, arrow)
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  # library(arrow)  # using fwrite instead
})

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NODE_PATH   <- file.path(BASE, "RNA-seq/results/network/network_nodes.csv")
COLOC_PATH  <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
OUTDIR      <- file.path(BASE, "RNA-seq/results/network")
OUT_CSV <- file.path(OUTDIR, "edges_genetic.csv")

PP4_THRESH  <- 0.5    # Minimum PP.H4.abf to count as colocalized
LOCUS_WINDOW <- 1e6   # 1 Mb window for single-linkage locus clustering

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 254: Genetic Co-Localization Edges ===\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n")
cat("PP.H4 threshold:", PP4_THRESH, "\n")
cat("Locus window:", format(LOCUS_WINDOW, big.mark = ","), "bp\n\n")

# ---------------------------------------------------------------------------
# 1. Load node set V
# ---------------------------------------------------------------------------
cat("--- Loading node set ---\n")
nodes <- fread(NODE_PATH)
cat("Node set V:", nrow(nodes), "genes\n")

# Build lookup: prefer ensembl_id, but also allow matching by gene symbol
node_ensembl <- unique(nodes$ensembl_id)
node_ensembl_unversioned <- unique(sub("\\.\\d+$", "", nodes$ensembl_id))
node_symbols <- unique(nodes$human_symbol[!is.na(nodes$human_symbol) & nodes$human_symbol != ""])
cat("Unique Ensembl IDs:", length(node_ensembl_unversioned), "\n")
cat("Unique gene symbols:", length(node_symbols), "\n\n")

# ---------------------------------------------------------------------------
# 2. Load aggregated COLOC results
# ---------------------------------------------------------------------------
cat("--- Loading SuSiE-COLOC results ---\n")
t_load <- Sys.time()

# T0.4 swap (review A13#6): prefer SuSiE PP4 over ABF (project-canonical; matches
# 207/208/210/213). Read PP.H4.susie when present and overwrite PP.H4.abf in place
# (SuSiE where available, ABF fallback) so the edge filter below uses the canonical PP4.
.have <- names(fread(COLOC_PATH, nrows = 0L))
.sel  <- intersect(c("gene", "ensembl", "chr", "gwas_name", "PP.H4.abf", "PP.H4.susie", "top_snp"), .have)
coloc <- fread(COLOC_PATH, select = .sel)
if ("PP.H4.susie" %in% names(coloc)) {
  coloc[, PP.H4.abf := fifelse(!is.na(PP.H4.susie), PP.H4.susie, PP.H4.abf)]
}
cat("Total gene-GWAS pairs:", format(nrow(coloc), big.mark = ","), "\n")
cat("GWAS studies:", uniqueN(coloc$gwas_name), "\n")
cat("Load time:", round(difftime(Sys.time(), t_load, units = "secs"), 1), "s\n\n")

# ---------------------------------------------------------------------------
# 3. Filter to significant colocalizations (PP.H4 > threshold)
# ---------------------------------------------------------------------------
cat("--- Filtering to PP.H4 >", PP4_THRESH, "---\n")
sig <- coloc[!is.na(PP.H4.abf) & PP.H4.abf > PP4_THRESH]
cat("Gene-GWAS pairs with PP.H4 >", PP4_THRESH, ":", format(nrow(sig), big.mark = ","), "\n")
cat("Unique genes:", uniqueN(sig$gene), "\n")
cat("GWAS with hits:", uniqueN(sig$gwas_name), "\n\n")

# ---------------------------------------------------------------------------
# 4. Parse top_snp position for locus clustering
# ---------------------------------------------------------------------------
cat("--- Parsing top_snp positions ---\n")

# top_snp format is "chr:position" (e.g., "1:109817590")
sig[, top_snp_pos := as.numeric(sub("^[^:]+:", "", top_snp))]
n_no_pos <- sum(is.na(sig$top_snp_pos))
if (n_no_pos > 0) {
  cat("WARNING:", n_no_pos, "entries with unparseable top_snp position, dropping\n")
  sig <- sig[!is.na(top_snp_pos)]
}
cat("Entries with valid position:", nrow(sig), "\n\n")

# ---------------------------------------------------------------------------
# 5. Assign locus IDs via single-linkage clustering within gwas+chr
# ---------------------------------------------------------------------------
# Within each GWAS + chromosome combination, cluster genes whose top_snp
# positions are within LOCUS_WINDOW of each other using single-linkage.
# ---------------------------------------------------------------------------
cat("--- Assigning locus IDs ---\n")
t_locus <- Sys.time()

sig[, gwas_chr := paste0(gwas_name, ":", chr)]

assign_loci <- function(positions, window) {
  # Single-linkage clustering: sort positions, merge neighbors within window
  n <- length(positions)
  if (n <= 1) return(rep(1L, n))

  ord <- order(positions)
  sorted <- positions[ord]
  cluster <- integer(n)
  cluster[ord[1]] <- 1L
  cid <- 1L

  for (i in 2:n) {
    if (sorted[i] - sorted[i - 1] > window) {
      cid <- cid + 1L
    }
    cluster[ord[i]] <- cid
  }
  return(cluster)
}

sig[, locus_local := assign_loci(top_snp_pos, LOCUS_WINDOW), by = gwas_chr]
sig[, locus_id := paste0(gwas_chr, ":L", locus_local)]

n_loci <- uniqueN(sig$locus_id)
cat("Total loci across all GWAS:", n_loci, "\n")

# Count genes per locus
locus_counts <- sig[, .(n_genes = uniqueN(gene)), by = locus_id]
n_multi <- sum(locus_counts$n_genes >= 2)
cat("Loci with 2+ genes:", n_multi, "\n")
cat("Loci with 3+ genes:", sum(locus_counts$n_genes >= 3), "\n")
cat("Max genes at a single locus:", max(locus_counts$n_genes), "\n")
cat("Locus assignment time:", round(difftime(Sys.time(), t_locus, units = "secs"), 1), "s\n\n")

# ---------------------------------------------------------------------------
# 6. Build edges: all pairs within each multi-gene locus
# ---------------------------------------------------------------------------
cat("--- Building pairwise edges ---\n")
t_edge <- Sys.time()

# Keep only loci with 2+ genes
multi_loci <- locus_counts[n_genes >= 2, locus_id]
sig_multi <- sig[locus_id %in% multi_loci]

cat("Gene-GWAS entries at multi-gene loci:", format(nrow(sig_multi), big.mark = ","), "\n")

# For each locus: create all gene pairs and compute edge weight
# Use data.table self-join within each locus
edge_list <- sig_multi[, {
  genes <- unique(gene)
  pp4s  <- PP.H4.abf[match(genes, gene)]
  n_g   <- length(genes)

  if (n_g < 2) {
    NULL
  } else {
    # All unique pairs (i < j by index)
    pairs <- combn(n_g, 2)
    data.table(
      gene_a    = genes[pairs[1, ]],
      gene_b    = genes[pairs[2, ]],
      pp4_a     = pp4s[pairs[1, ]],
      pp4_b     = pp4s[pairs[2, ]],
      shared_gwas = gwas_name[1]
    )
  }
}, by = locus_id]

cat("Raw edges (before dedup):", format(nrow(edge_list), big.mark = ","), "\n")

# Edge weight = min(PP.H4_a, PP.H4_b) -- bottleneck confidence
edge_list[, raw_score := pmin(pp4_a, pp4_b)]

cat("Edge build time:", round(difftime(Sys.time(), t_edge, units = "secs"), 1), "s\n\n")

# ---------------------------------------------------------------------------
# 7. Ensure canonical ordering: gene_a < gene_b alphabetically
# ---------------------------------------------------------------------------
cat("--- Canonical ordering and dedup ---\n")

swap <- edge_list$gene_a > edge_list$gene_b
if (any(swap)) {
  tmp_a <- edge_list$gene_a[swap]
  tmp_pp4_a <- edge_list$pp4_a[swap]
  edge_list[swap, gene_a := gene_b]
  edge_list[swap, pp4_a := pp4_b]
  edge_list[swap, gene_b := tmp_a]
  edge_list[swap, pp4_b := tmp_pp4_a]
  rm(tmp_a, tmp_pp4_a)
}

# Remove self-loops (shouldn't exist but safety check)
edge_list <- edge_list[gene_a != gene_b]

# If a gene pair shares multiple loci (same or different GWAS), keep the MAX weight
# Concatenate all shared_gwas and locus_id for metadata
edges <- edge_list[, .(
  raw_score  = max(raw_score),
  shared_gwas = paste(unique(shared_gwas), collapse = ";"),
  locus_id    = paste(unique(locus_id), collapse = ";"),
  n_shared_loci = uniqueN(locus_id)
), by = .(gene_a, gene_b)]

cat("Edges after max-dedup:", format(nrow(edges), big.mark = ","), "\n")
cat("Pairs with 2+ shared loci:", sum(edges$n_shared_loci >= 2), "\n\n")

# ---------------------------------------------------------------------------
# 8. Filter to genes in node set V
# ---------------------------------------------------------------------------
cat("--- Filtering to node set V ---\n")

# Match by ensembl_id (unversioned) or gene symbol
# The COLOC file uses gene symbols in the 'gene' column
in_V <- function(g) g %in% node_symbols | g %in% node_ensembl_unversioned

mask <- in_V(edges$gene_a) & in_V(edges$gene_b)
cat("Edges before V filter:", format(nrow(edges), big.mark = ","), "\n")
edges <- edges[mask]
cat("Edges after V filter:", format(nrow(edges), big.mark = ","), "\n\n")

# ---------------------------------------------------------------------------
# 9. Format metadata and save as parquet
# ---------------------------------------------------------------------------
cat("--- Saving parquet ---\n")

edges[, metadata := paste0("shared_gwas=", shared_gwas, ";locus_id=", locus_id,
                           ";n_shared_loci=", n_shared_loci)]

out <- edges[, .(gene_a, gene_b, raw_score, metadata)]
setorder(out, -raw_score)

fwrite(out, OUT_CSV)
cat("Saved:", OUT_CSV, "\n")
cat("File size:", round(file.size(OUT_CSV) / 1e6, 2), "MB\n\n")

# ---------------------------------------------------------------------------
# 10. Summary statistics
# ---------------------------------------------------------------------------
cat("=== Summary ===\n")
cat("Loci with 2+ colocalized genes:", n_multi, "\n")
cat("Total edges:", format(nrow(out), big.mark = ","), "\n")
cat("Mean weight (min PP.H4):", round(mean(out$raw_score), 4), "\n")
cat("Median weight:", round(median(out$raw_score), 4), "\n")
cat("Min weight:", round(min(out$raw_score), 4), "\n")
cat("Max weight:", round(max(out$raw_score), 4), "\n")

# Weight distribution
cat("\n--- Weight distribution ---\n")
breaks <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0)
for (i in seq_along(breaks[-1])) {
  lo <- breaks[i]
  hi <- breaks[i + 1]
  n_in <- sum(out$raw_score >= lo & out$raw_score < hi)
  cat(sprintf("  [%.2f, %.2f): %d edges\n", lo, hi, n_in))
}
cat(sprintf("  = 1.0:       %d edges\n", sum(out$raw_score == 1.0)))

# Top 10 highest-degree genes
cat("\n--- Top 10 highest-degree genes ---\n")
degree_a <- out[, .N, by = gene_a]
degree_b <- out[, .N, by = gene_b]
setnames(degree_a, c("gene", "deg_a"))
setnames(degree_b, c("gene", "deg_b"))
degree <- merge(degree_a, degree_b, by = "gene", all = TRUE)
degree[is.na(deg_a), deg_a := 0L]
degree[is.na(deg_b), deg_b := 0L]
degree[, degree := deg_a + deg_b]
setorder(degree, -degree)

cat(sprintf("%-18s %s\n", "gene", "degree"))
cat(strrep("-", 30), "\n")
for (i in seq_len(min(10, nrow(degree)))) {
  cat(sprintf("%-18s %d\n", degree$gene[i], degree$degree[i]))
}

# Top 10 highest-weight edges
cat("\n--- Top 10 highest-weight edges ---\n")
cat(sprintf("%-15s %-15s %8s  %s\n", "gene_a", "gene_b", "weight", "gwas"))
cat(strrep("-", 60), "\n")
for (i in seq_len(min(10, nrow(out)))) {
  gwas_short <- sub(";locus_id=.*", "", sub("^shared_gwas=", "", out$metadata[i]))
  cat(sprintf("%-15s %-15s %8.4f  %s\n",
              out$gene_a[i], out$gene_b[i], out$raw_score[i], gwas_short))
}

t1 <- Sys.time()
cat("\nElapsed:", round(difftime(t1, t0, units = "secs"), 1), "s\n")
cat("Done.\n")
