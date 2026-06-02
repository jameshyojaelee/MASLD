#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_251_ppi
#SBATCH --output=logs/net_251_ppi_%j.out
#SBATCH --error=logs/net_251_ppi_%j.err
# ===========================================================================
# Script 251: Build PPI edges from STRING v12
# ===========================================================================
# Purpose: Extract protein-protein interaction edges from STRING v12 for the
#          Bayesian multiplex network. Filters to the node set defined by
#          Script 250 (network_nodes.csv).
#
# Input:
#   - RNA-seq/results/network/network_nodes.csv  (node set, human_symbol col)
#   - data/string_ppi/9606.protein.links.v12.0.txt.gz  (STRING v12 links)
#   - data/string_ppi/9606.protein.info.v12.0.txt.gz   (ENSP -> symbol map)
#
# Output:
#   - RNA-seq/results/network/edges_ppi.csv
#     Columns: gene_a, gene_b, raw_score, metadata
#
# Logic:
#   1. Load network node set (V) from Script 250
#   2. Load STRING v12, filter combined_score > 400
#   3. Map ENSP IDs to gene symbols via STRING info file
#   4. Filter to gene pairs where BOTH genes are in V
#   5. Remove self-loops, deduplicate (max score per pair), gene_a < gene_b
#   6. Save as parquet
#
# Environment: rnaseq (R 4.4+, data.table, arrow)
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  # library(arrow)  # using fwrite instead
})

t0 <- Sys.time()
cat("=== Script 251: PPI Edges from STRING v12 ===\n")
cat("Start:", format(t0), "\n")

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

node_file   <- file.path(BASE, "RNA-seq/results/network/network_nodes.csv")
ppi_file    <- file.path(BASE, "data/string_ppi/9606.protein.links.v12.0.txt.gz")
info_file   <- file.path(BASE, "data/string_ppi/9606.protein.info.v12.0.txt.gz")
out_csv <- file.path(BASE, "RNA-seq/results/network/edges_ppi.csv")

stopifnot(file.exists(node_file))
stopifnot(file.exists(ppi_file))
stopifnot(file.exists(info_file))

# ---------------------------------------------------------------------------
# 1. Load network node set
# ---------------------------------------------------------------------------
cat("\n--- Loading network nodes ---\n")
nodes <- fread(node_file)
cat("Node file:", nrow(nodes), "rows x", ncol(nodes), "columns\n")
cat("Columns:", paste(names(nodes), collapse = ", "), "\n")

V <- unique(nodes$human_symbol)
V <- V[!is.na(V) & V != ""]
cat("Unique gene symbols in V:", length(V), "\n")

# ---------------------------------------------------------------------------
# 2. Load STRING protein info (ENSP -> symbol mapping)
# ---------------------------------------------------------------------------
cat("\n--- Loading STRING protein info ---\n")
ppi_info <- fread(info_file)
cat("STRING info rows:", nrow(ppi_info), "\n")

# Strip "9606." prefix from STRING protein IDs (same approach as Script 59)
ppi_info[, ensp := sub("^9606\\.", "", `#string_protein_id`)]
ensp2sym <- setNames(ppi_info$preferred_name, ppi_info$ensp)
cat("ENSP-to-symbol mappings:", length(ensp2sym), "\n")

# ---------------------------------------------------------------------------
# 3. Load STRING links, filter combined_score > 400
# ---------------------------------------------------------------------------
cat("\n--- Loading STRING links ---\n")
ppi <- fread(ppi_file)
cat("Raw STRING edges:", format(nrow(ppi), big.mark = ","), "\n")

ppi <- ppi[combined_score > 400]
cat("Edges after score > 400 filter:", format(nrow(ppi), big.mark = ","), "\n")

# ---------------------------------------------------------------------------
# 4. Map ENSP IDs to gene symbols
# ---------------------------------------------------------------------------
cat("\n--- Mapping ENSP to gene symbols ---\n")
ppi[, p1 := sub("^9606\\.", "", protein1)]
ppi[, p2 := sub("^9606\\.", "", protein2)]
ppi[, sym1 := ensp2sym[p1]]
ppi[, sym2 := ensp2sym[p2]]

# Drop edges where either protein has no symbol mapping
n_before <- nrow(ppi)
ppi <- ppi[!is.na(sym1) & !is.na(sym2)]
cat("Edges with both symbols mapped:", format(nrow(ppi), big.mark = ","),
    sprintf(" (dropped %s unmapped)\n", format(n_before - nrow(ppi), big.mark = ",")))

# ---------------------------------------------------------------------------
# 5. Filter to gene pairs where BOTH genes are in V
# ---------------------------------------------------------------------------
cat("\n--- Filtering to node set V ---\n")
ppi <- ppi[sym1 %chin% V & sym2 %chin% V]
cat("Edges with both genes in V:", format(nrow(ppi), big.mark = ","), "\n")

# ---------------------------------------------------------------------------
# 6. Remove self-loops, deduplicate, ensure gene_a < gene_b
# ---------------------------------------------------------------------------
cat("\n--- Deduplicating ---\n")

# Remove self-loops
ppi <- ppi[sym1 != sym2]
cat("After removing self-loops:", format(nrow(ppi), big.mark = ","), "\n")

# Normalize pair order: gene_a < gene_b alphabetically
ppi[, gene_a := fifelse(sym1 < sym2, sym1, sym2)]
ppi[, gene_b := fifelse(sym1 < sym2, sym2, sym1)]

# Deduplicate: keep max combined_score per pair
ppi_dedup <- ppi[, .(combined_score = max(combined_score)), by = .(gene_a, gene_b)]
cat("After deduplication (max score per pair):", format(nrow(ppi_dedup), big.mark = ","), "\n")

# ---------------------------------------------------------------------------
# 7. Build output table
# ---------------------------------------------------------------------------
cat("\n--- Building output ---\n")

edges <- ppi_dedup[, .(
  gene_a    = gene_a,
  gene_b    = gene_b,
  raw_score = combined_score / 1000,
  metadata  = "combined_score_only"
)]

# ---------------------------------------------------------------------------
# 8. Save as parquet
# ---------------------------------------------------------------------------
cat("\n--- Saving parquet ---\n")
fwrite(edges, out_csv)
cat("Saved:", out_csv, "\n")
cat("File size:", sprintf("%.1f MB", file.size(out_csv) / 1e6), "\n")

# ---------------------------------------------------------------------------
# 9. Summary statistics
# ---------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Total edges:", format(nrow(edges), big.mark = ","), "\n")
cat("Mean raw_score:", round(mean(edges$raw_score), 3), "\n")
cat("Median raw_score:", round(median(edges$raw_score), 3), "\n")
cat("Score range:", round(min(edges$raw_score), 3), "-",
    round(max(edges$raw_score), 3), "\n")

# Top 10 highest-degree genes
degree <- rbind(
  edges[, .(gene = gene_a)],
  edges[, .(gene = gene_b)]
)[, .N, by = gene][order(-N)]

cat("\nTop 10 highest-degree genes:\n")
print(degree[1:min(10, nrow(degree))])

cat("\nDegree distribution quantiles:\n")
cat(sprintf("  Min: %d  Q25: %d  Median: %d  Q75: %d  Max: %d\n",
            min(degree$N), quantile(degree$N, 0.25),
            median(degree$N), quantile(degree$N, 0.75),
            max(degree$N)))

cat("\nUnique genes in edges:", nrow(degree), "\n")

t1 <- Sys.time()
cat("\nElapsed:", format(round(difftime(t1, t0, units = "mins"), 1)), "\n")
cat("Done:", format(t1), "\n")
