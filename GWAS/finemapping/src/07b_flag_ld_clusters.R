#!/usr/bin/env Rscript
# 07b_flag_ld_clusters.R
# Identifies LD contamination clusters in COLOC results where multiple adjacent
# genes share the same or nearby lead SNP and all receive high PP.H4.
#
# These clusters represent a single locus signal that should not be counted as
# independent gene-level evidence.
#
# Usage: Rscript 07b_flag_ld_clusters.R
# Environment: finemapping conda env

library(data.table)
library(dplyr)

# --- Configuration -----------------------------------------------------------

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
COLOC_DIR <- file.path(FM_DIR, "results/susie_coloc")

# Clustering parameters
PP4_THRESHOLD     <- 0.5   # Only consider genes with PP.H4 above this
CLUSTER_WINDOW_KB <- 100   # Group genes whose top_snp positions are within this window (kb)
SHARED_SNP_KB     <- 10    # Flag clusters where top_snps are within this distance as "shared SNP"
MIN_CLUSTER_SIZE  <- 3     # Minimum genes in a cluster to flag

# Known biology: for each locus, the most likely causal gene based on
# functional evidence (eQTL fine-mapping, CRISPR, known biology).
# Keys are "chr:approx_mb" to allow fuzzy matching.
KNOWN_CAUSAL <- list(
  "1:109" = "SORT1",       # SORT1 locus: Musunuru et al. 2010, CRISPR-validated
  "6:31"  = "VARS2",       # HLA region: no single causal gene, use highest PP4
  "6:32"  = "VARS2",       # HLA extended
  "8:9"   = "PPP1R3B",     # PPP1R3B: functional glycogen phosphatase
  "2:27"  = "C2orf16"      # C2orf16/GPN1/ATRAID locus
)

# --- Input files -------------------------------------------------------------

all_gwas_file <- file.path(COLOC_DIR, "susie_coloc_all_gwas.csv")
gene_level_file <- file.path(COLOC_DIR, "gene_level_coloc.csv")

if (!file.exists(all_gwas_file)) {
  stop("Missing: ", all_gwas_file, "\nRun 07_combine_susie_coloc.R first.")
}
if (!file.exists(gene_level_file)) {
  stop("Missing: ", gene_level_file, "\nRun 07_combine_susie_coloc.R first.")
}

cat("============================================================\n")
cat("LD contamination cluster detection\n")
cat("============================================================\n")
cat("Parameters:\n")
cat("  PP.H4 threshold:     ", PP4_THRESHOLD, "\n")
cat("  Cluster window:      ", CLUSTER_WINDOW_KB, "kb\n")
cat("  Shared SNP window:   ", SHARED_SNP_KB, "kb\n")
cat("  Min cluster size:    ", MIN_CLUSTER_SIZE, "\n")
cat("============================================================\n\n")

# --- Load data ---------------------------------------------------------------

all_gwas <- fread(all_gwas_file)
cat("Loaded", nrow(all_gwas), "gene-GWAS entries from", all_gwas_file, "\n")

# Filter to significant colocalisations
sig <- all_gwas[PP.H4.abf >= PP4_THRESHOLD & !is.na(top_snp) & top_snp != ""]
cat("Entries with PP.H4 >=", PP4_THRESHOLD, ":", nrow(sig), "\n\n")

if (nrow(sig) == 0) {
  cat("No significant colocalisations found. Exiting.\n")
  quit(save = "no", status = 0)
}

# Parse top_snp into chr and position
sig[, c("snp_chr", "snp_pos") := {
  parts <- tstrsplit(top_snp, ":", fixed = TRUE)
  list(as.integer(parts[[1]]), as.integer(parts[[2]]))
}]

# Drop rows where parsing failed
sig <- sig[!is.na(snp_chr) & !is.na(snp_pos)]
cat("After parsing top_snp:", nrow(sig), "entries\n\n")

# --- Cluster detection -------------------------------------------------------
# For each GWAS separately, group genes on the same chromosome whose
# top_snp positions fall within CLUSTER_WINDOW_KB of each other.
# Uses single-linkage clustering via hierarchical clustering.

cluster_results <- list()
cluster_counter <- 0

gwas_names <- unique(sig$gwas_name)
cat("Processing", length(gwas_names), "GWAS...\n")

for (gw in gwas_names) {
  gw_data <- sig[gwas_name == gw]
  chromosomes <- unique(gw_data$snp_chr)

  for (chr_val in chromosomes) {
    chr_data <- gw_data[snp_chr == chr_val]

    if (nrow(chr_data) < MIN_CLUSTER_SIZE) next

    # Single-linkage hierarchical clustering on SNP positions
    positions <- chr_data$snp_pos
    if (length(positions) < 2) next

    dist_mat <- dist(positions)
    hc <- hclust(dist_mat, method = "single")
    chr_data[, cluster := cutree(hc, h = CLUSTER_WINDOW_KB * 1000)]

    # Process each cluster
    for (cl in unique(chr_data$cluster)) {
      cl_data <- chr_data[cluster == cl]

      if (nrow(cl_data) < MIN_CLUSTER_SIZE) next

      cluster_counter <- cluster_counter + 1

      # Determine position range
      pos_min <- min(cl_data$snp_pos)
      pos_max <- max(cl_data$snp_pos)
      pos_range_kb <- (pos_max - pos_min) / 1000

      # Check if SNPs are truly shared (within SHARED_SNP_KB)
      # Compute pairwise max distance between all top_snps in the cluster
      unique_positions <- sort(unique(cl_data$snp_pos))
      snp_spread_kb <- (max(unique_positions) - min(unique_positions)) / 1000
      shared_snp_flag <- snp_spread_kb <= SHARED_SNP_KB

      # Determine the "shared" SNP — the mode (most frequent top_snp)
      snp_freq <- sort(table(cl_data$top_snp), decreasing = TRUE)
      dominant_snp <- names(snp_freq)[1]
      n_with_dominant <- as.integer(snp_freq[1])

      # Identify best gene by PP.H4
      best_idx <- which.max(cl_data$PP.H4.abf)
      best_gene_pp4 <- cl_data$gene[best_idx]
      best_pp4 <- cl_data$PP.H4.abf[best_idx]

      # Check known biology for override
      chr_mb <- paste0(chr_val, ":", floor(pos_min / 1e6))
      known_gene <- KNOWN_CAUSAL[[chr_mb]]
      if (!is.null(known_gene) && known_gene %in% cl_data$gene) {
        best_gene <- known_gene
        selection_method <- "known_biology"
      } else {
        best_gene <- best_gene_pp4
        selection_method <- "highest_pp4"
      }

      # Classify the flag reason
      if (shared_snp_flag && nrow(cl_data) >= 5) {
        flag_reason <- "HLA-like_shared_snp_large_cluster"
      } else if (shared_snp_flag) {
        flag_reason <- "shared_snp_cluster"
      } else if (pos_range_kb <= 500) {
        flag_reason <- "ld_proximity_cluster"
      } else {
        flag_reason <- "extended_ld_region"
      }

      # Gene list sorted by PP.H4 descending
      gene_order <- cl_data[order(-PP.H4.abf)]
      genes_str <- paste(gene_order$gene, collapse = ", ")

      cluster_results[[cluster_counter]] <- data.table(
        cluster_id     = sprintf("LD_%03d", cluster_counter),
        chr            = chr_val,
        position_range = sprintf("%s:%d-%d (%.0fkb)",
                                 chr_val, pos_min, pos_max, pos_range_kb),
        n_genes        = nrow(cl_data),
        genes          = genes_str,
        shared_snp     = dominant_snp,
        n_sharing_top_snp = n_with_dominant,
        snp_spread_kb  = round(snp_spread_kb, 2),
        best_gene      = best_gene,
        best_pp4       = round(best_pp4, 4),
        selection_method = selection_method,
        gwas_name      = gw,
        flag_reason    = flag_reason
      )
    }
  }
}

# --- Assemble and write output -----------------------------------------------

if (length(cluster_results) == 0) {
  cat("\nNo LD contamination clusters found.\n")
  # Write empty file with header
  empty_dt <- data.table(
    cluster_id = character(), chr = integer(),
    position_range = character(), n_genes = integer(),
    genes = character(), shared_snp = character(),
    n_sharing_top_snp = integer(), snp_spread_kb = numeric(),
    best_gene = character(), best_pp4 = numeric(),
    selection_method = character(), gwas_name = character(),
    flag_reason = character()
  )
  out_file <- file.path(COLOC_DIR, "ld_contamination_clusters.csv")
  fwrite(empty_dt, out_file)
  cat("Written empty table to:", out_file, "\n")
  quit(save = "no", status = 0)
}

clusters <- rbindlist(cluster_results)
clusters <- clusters[order(chr, as.integer(gsub(".*:(\\d+)-.*", "\\1", position_range)), gwas_name)]

out_file <- file.path(COLOC_DIR, "ld_contamination_clusters.csv")
fwrite(clusters, out_file)

# --- Summary -----------------------------------------------------------------

cat("\n============================================================\n")
cat("LD contamination cluster summary\n")
cat("============================================================\n")
cat("Total clusters found:", nrow(clusters), "\n")
cat("Unique loci (chr:range):", length(unique(clusters$position_range)), "\n")
cat("GWAS with clusters:", length(unique(clusters$gwas_name)), "\n\n")

cat("--- By flag reason ---\n")
reason_summary <- clusters[, .(n_clusters = .N, total_genes = sum(n_genes)),
                           by = flag_reason][order(-n_clusters)]
for (i in seq_len(nrow(reason_summary))) {
  cat(sprintf("  %-40s %3d clusters, %3d genes\n",
              reason_summary$flag_reason[i],
              reason_summary$n_clusters[i],
              reason_summary$total_genes[i]))
}

cat("\n--- Largest clusters ---\n")
top_clusters <- clusters[order(-n_genes)][1:min(10, nrow(clusters))]
for (i in seq_len(nrow(top_clusters))) {
  row <- top_clusters[i]
  cat(sprintf("  %s | %s | %d genes | best=%s (PP4=%.3f) | %s\n",
              row$cluster_id, row$position_range, row$n_genes,
              row$best_gene, row$best_pp4, row$gwas_name))
  cat(sprintf("    genes: %s\n", row$genes))
}

# --- Validation: check known clusters ----------------------------------------

cat("\n--- Validation: known clusters ---\n")

# HLA chr6 ~31.7Mb
hla_hits <- clusters[chr == 6 & grepl("317[0-9]{5}", position_range)]
if (nrow(hla_hits) > 0) {
  cat("  [PASS] HLA cluster found:", nrow(hla_hits), "GWAS\n")
  for (i in seq_len(nrow(hla_hits))) {
    cat(sprintf("    %s: %d genes (%s)\n",
                hla_hits$gwas_name[i], hla_hits$n_genes[i], hla_hits$genes[i]))
  }
} else {
  cat("  [WARN] HLA cluster NOT found\n")
}

# SORT1 chr1 ~109.8Mb
sort1_hits <- clusters[chr == 1 & grepl("109[0-9]{6}", position_range)]
if (nrow(sort1_hits) > 0) {
  cat("  [PASS] SORT1 cluster found:", nrow(sort1_hits), "GWAS\n")
  for (i in seq_len(nrow(sort1_hits))) {
    cat(sprintf("    %s: %d genes (%s)\n",
                sort1_hits$gwas_name[i], sort1_hits$n_genes[i], sort1_hits$genes[i]))
  }
} else {
  cat("  [WARN] SORT1 cluster NOT found\n")
}

# PPP1R3B chr8 ~9.2Mb — only 2 genes, below MIN_CLUSTER_SIZE=3, so note it
ppp_hits <- clusters[chr == 8 & grepl("91[0-9]{5}", position_range)]
if (nrow(ppp_hits) > 0) {
  cat("  [PASS] PPP1R3B cluster found:", nrow(ppp_hits), "GWAS\n")
} else {
  cat("  [NOTE] PPP1R3B cluster not found (expected: only 2 genes, below min_cluster_size=3)\n")
}

# C2orf16/GPN1 chr2 ~27Mb
c2_hits <- clusters[chr == 2 & grepl("27[0-9]{6}", position_range)]
if (nrow(c2_hits) > 0) {
  cat("  [PASS] C2orf16/GPN1 cluster found:", nrow(c2_hits), "GWAS\n")
  for (i in seq_len(nrow(c2_hits))) {
    cat(sprintf("    %s: %d genes (%s)\n",
                c2_hits$gwas_name[i], c2_hits$n_genes[i], c2_hits$genes[i]))
  }
} else {
  cat("  [WARN] C2orf16/GPN1 cluster NOT found\n")
}

# --- Summary of total gene-GWAS entries flagged ------------------------------

total_flagged <- sum(clusters$n_genes)
total_sig <- nrow(sig)
cat(sprintf("\n--- Impact ---\n"))
cat(sprintf("  Total gene-GWAS entries with PP.H4 >= %.1f: %d\n", PP4_THRESHOLD, total_sig))
cat(sprintf("  Gene-GWAS entries in LD clusters: %d (%.1f%%)\n",
            total_flagged, 100 * total_flagged / total_sig))
cat(sprintf("  After deduplication (keeping best per cluster): %d entries removed\n",
            total_flagged - nrow(clusters)))

cat("\nWritten to:", out_file, "\n")
cat("============================================================\n")
