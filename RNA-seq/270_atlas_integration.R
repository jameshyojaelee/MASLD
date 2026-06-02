#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_270_atlas
#SBATCH --output=logs/net_270_atlas_%j.out
#SBATCH --error=logs/net_270_atlas_%j.err
# ===========================================================================
# Script 270: Atlas Integration — Network-Derived Columns
# ===========================================================================
# Purpose: Stage 4b of the Bayesian multiplex gene network pipeline. Adds 6
#          network-derived columns to the multi-evidence atlas.
#
# New columns:
#   network_community_macro       — macro community ID
#   network_community_meso        — meso community ID
#   network_community_micro       — micro community ID
#   network_degree_composite      — number of composite edges (P > 0.5)
#   network_top_neighbor          — gene symbol of highest-weighted neighbor
#   network_edge_multiplicity_max — max K_ij across all neighbors
#
# Input:
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#   - RNA-seq/results/network/communities/community_assignments.csv
#   - RNA-seq/results/network/composite_edges.parquet
#
# Output:
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (updated in-place)
#
# Follows the idempotent re-run guard pattern from Scripts 75 and 217:
# existing network_ columns are dropped before re-merge.
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

cat("=== Script 270: Network Atlas Integration ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS_FILE <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
NET_DIR    <- file.path(BASE, "RNA-seq/results/network")
COMM_FILE  <- file.path(NET_DIR, "communities/community_assignments.csv")
EDGES_FILE <- file.path(NET_DIR, "composite_edges.csv")

stopifnot(file.exists(ATLAS_FILE))

# ===========================================================================
# 1. Load atlas
# ===========================================================================
atlas <- fread(ATLAS_FILE)
n_col_start <- ncol(atlas)
cat(sprintf("Loaded atlas: %d genes x %d columns\n", nrow(atlas), n_col_start))

# Identify gene column
gene_col <- "human_symbol"
if (!gene_col %in% names(atlas)) {
  gene_col <- intersect(names(atlas), c("gene", "gene_symbol", "symbol"))[1]
  if (is.na(gene_col)) stop("No gene symbol column found in atlas")
}

# ===========================================================================
# 2. Re-run guard: remove existing network_ columns (idempotency)
# ===========================================================================
network_cols <- c(
  "network_community_macro",
  "network_community_meso",
  "network_community_micro",
  "network_degree_composite",
  "network_top_neighbor",
  "network_edge_multiplicity_max"
)

existing <- intersect(network_cols, names(atlas))
if (length(existing) > 0L) {
  cat(sprintf("Re-run detected: removing %d existing network columns: %s\n",
              length(existing), paste(existing, collapse = ", ")))
  atlas[, (existing) := NULL]
}

# ===========================================================================
# 3. Load community assignments
# ===========================================================================
cat("\n--- Loading community assignments ---\n")

if (file.exists(COMM_FILE)) {
  comm <- fread(COMM_FILE)
  cat(sprintf("  Community assignments: %d genes x %d columns\n", nrow(comm), ncol(comm)))

  # Identify gene column in community file
  comm_gene_col <- intersect(names(comm), c("gene", "human_symbol", "symbol",
                                             "node", "gene_symbol"))[1]
  if (is.na(comm_gene_col)) {
    cat("  WARNING: No gene column found in community_assignments.csv\n")
    cat(sprintf("    Available columns: %s\n", paste(head(names(comm), 10), collapse = ", ")))
    comm <- NULL
  } else {
    if (comm_gene_col != gene_col) {
      setnames(comm, comm_gene_col, gene_col)
    }

    # Identify community level columns (flexible naming)
    macro_col <- intersect(names(comm), c("community_macro", "macro", "macro_id"))[1]
    meso_col  <- intersect(names(comm), c("community_meso", "meso", "meso_id"))[1]
    micro_col <- intersect(names(comm), c("community_micro", "micro", "micro_id"))[1]

    # Build standardized community table
    comm_dt <- comm[, .SD[1], by = c(gene_col)]  # deduplicate
    cols_to_keep <- gene_col

    if (!is.na(macro_col)) {
      comm_dt[, network_community_macro := as.integer(get(macro_col))]
      cols_to_keep <- c(cols_to_keep, "network_community_macro")
    }
    if (!is.na(meso_col)) {
      comm_dt[, network_community_meso := as.integer(get(meso_col))]
      cols_to_keep <- c(cols_to_keep, "network_community_meso")
    }
    if (!is.na(micro_col)) {
      comm_dt[, network_community_micro := as.integer(get(micro_col))]
      cols_to_keep <- c(cols_to_keep, "network_community_micro")
    }

    comm_dt <- comm_dt[, ..cols_to_keep]
    cat(sprintf("  Prepared community columns: %s\n",
                paste(setdiff(cols_to_keep, gene_col), collapse = ", ")))
  }
} else {
  cat(sprintf("  WARNING: %s not found — community columns will be NA\n", COMM_FILE))
  comm <- NULL
  comm_dt <- NULL
}

# ===========================================================================
# 4. Load composite edges and compute per-gene network metrics
# ===========================================================================
cat("\n--- Loading composite edges ---\n")

if (file.exists(EDGES_FILE)) {
  edges <- fread(EDGES_FILE)
  cat(sprintf("  Composite edges: %d rows x %d columns\n", nrow(edges), ncol(edges)))

  # Identify columns
  ga_col <- intersect(names(edges), c("gene_a", "source", "node_a", "from"))[1]
  gb_col <- intersect(names(edges), c("gene_b", "target", "node_b", "to"))[1]
  # Posterior probability or composite weight
  prob_col <- intersect(names(edges), c("convergence_score", "posterior_prob", "prob", "weight",
                                         "composite_weight", "P", "p"))[1]
  # Layer count / multiplicity
  k_col <- intersect(names(edges), c("K_ij", "k_ij", "n_layers", "multiplicity",
                                      "layer_count", "edge_multiplicity"))[1]

  if (is.na(ga_col) || is.na(gb_col)) {
    cat(sprintf("  WARNING: Cannot identify gene_a/gene_b columns (found: %s)\n",
                paste(head(names(edges), 10), collapse = ", ")))
    edge_metrics <- NULL
  } else {
    cat(sprintf("  Using columns: gene_a=%s, gene_b=%s", ga_col, gb_col))
    if (!is.na(prob_col)) cat(sprintf(", prob=%s", prob_col))
    if (!is.na(k_col))    cat(sprintf(", K=%s", k_col))
    cat("\n")

    # Filter to P > 0.5 if probability column exists
    if (!is.na(prob_col)) {
      n_before <- nrow(edges)
      edges <- edges[get(prob_col) > 0.5]
      cat(sprintf("  Filtered to P > 0.5: %d -> %d edges\n", n_before, nrow(edges)))
    }

    if (nrow(edges) == 0) {
      cat("  WARNING: No edges pass P > 0.5 threshold\n")
      edge_metrics <- NULL
    } else {
      # --- Degree: count edges per gene ---
      deg_a <- edges[, .(n = .N), by = c(ga_col)]
      setnames(deg_a, ga_col, "gene_tmp")
      deg_b <- edges[, .(n = .N), by = c(gb_col)]
      setnames(deg_b, gb_col, "gene_tmp")
      degree <- rbind(deg_a, deg_b)[, .(network_degree_composite = sum(n)), by = gene_tmp]
      setnames(degree, "gene_tmp", gene_col)

      # --- Top neighbor: highest weight per gene ---
      # For gene_a side
      if (!is.na(prob_col)) {
        weight_col_name <- prob_col
      } else {
        # Use constant weight if no weight column
        edges[, const_weight := 1.0]
        weight_col_name <- "const_weight"
      }

      top_a <- edges[, .SD[which.max(get(weight_col_name))], by = c(ga_col)]
      top_a <- top_a[, .(gene_tmp = get(ga_col), neighbor = get(gb_col),
                          weight = get(weight_col_name))]
      top_b <- edges[, .SD[which.max(get(weight_col_name))], by = c(gb_col)]
      top_b <- top_b[, .(gene_tmp = get(gb_col), neighbor = get(ga_col),
                          weight = get(weight_col_name))]
      top_both <- rbind(top_a, top_b)
      top_neighbor <- top_both[, .SD[which.max(weight)], by = gene_tmp]
      top_neighbor <- top_neighbor[, .(gene_tmp, network_top_neighbor = neighbor)]
      setnames(top_neighbor, "gene_tmp", gene_col)

      # --- Max multiplicity: max K_ij across all neighbors ---
      if (!is.na(k_col)) {
        kmax_a <- edges[, .(k = max(get(k_col), na.rm = TRUE)), by = c(ga_col)]
        setnames(kmax_a, ga_col, "gene_tmp")
        kmax_b <- edges[, .(k = max(get(k_col), na.rm = TRUE)), by = c(gb_col)]
        setnames(kmax_b, gb_col, "gene_tmp")
        kmax <- rbind(kmax_a, kmax_b)[, .(network_edge_multiplicity_max = max(k)), by = gene_tmp]
        setnames(kmax, "gene_tmp", gene_col)
      } else {
        # If no multiplicity column, set to NA
        kmax <- degree[, .(network_edge_multiplicity_max = NA_integer_), by = c(gene_col)]
      }

      # Merge edge metrics
      edge_metrics <- merge(degree, top_neighbor, by = gene_col, all = TRUE)
      edge_metrics <- merge(edge_metrics, kmax, by = gene_col, all = TRUE)
      cat(sprintf("  Edge metrics computed for %d genes\n", nrow(edge_metrics)))
    }
  }
} else {
  cat(sprintf("  WARNING: %s not found — edge metrics will be NA\n", EDGES_FILE))
  edge_metrics <- NULL
}

# ===========================================================================
# 5. Merge into atlas
# ===========================================================================
cat("\n--- Merging into atlas ---\n")

if (!is.null(comm_dt)) {
  atlas <- merge(atlas, comm_dt, by = gene_col, all.x = TRUE)
  cat(sprintf("  Merged community assignments: %d genes with values\n",
              sum(!is.na(atlas$network_community_macro))))
}

if (!is.null(edge_metrics)) {
  atlas <- merge(atlas, edge_metrics, by = gene_col, all.x = TRUE)
  cat(sprintf("  Merged edge metrics: %d genes with degree data\n",
              sum(!is.na(atlas$network_degree_composite))))
}

# Ensure all 6 columns exist (initialize as NA if source data was missing)
for (col in network_cols) {
  if (!col %in% names(atlas)) {
    cat(sprintf("  Initializing missing column %s as NA\n", col))
    if (grepl("community|multiplicity|degree", col)) {
      atlas[, (col) := NA_integer_]
    } else {
      atlas[, (col) := NA_character_]
    }
  }
}

# ===========================================================================
# 6. Save updated atlas
# ===========================================================================
cat("\n--- Writing updated atlas ---\n")

# Backup current atlas
backup_fp <- sub("\\.csv$",
                 sprintf("_backup_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
                 ATLAS_FILE)
file.copy(ATLAS_FILE, backup_fp)
cat(sprintf("  Backup: %s\n", backup_fp))

fwrite(atlas, ATLAS_FILE)
n_col_end <- ncol(atlas)
cat(sprintf("  Updated atlas: %d genes x %d columns (%d -> %d, +%d)\n",
            nrow(atlas), n_col_end, n_col_start, n_col_end, n_col_end - n_col_start))

# ===========================================================================
# 7. Summary statistics
# ===========================================================================
cat("\n--- Summary Statistics ---\n")

# Column coverage
for (col in network_cols) {
  n_nona <- sum(!is.na(atlas[[col]]))
  cat(sprintf("  %-40s: %d non-NA (%.1f%%)\n", col, n_nona,
              100 * n_nona / nrow(atlas)))
}

# Degree distribution
if (sum(!is.na(atlas$network_degree_composite)) > 0) {
  deg <- atlas$network_degree_composite[!is.na(atlas$network_degree_composite)]
  cat(sprintf("\n  Degree distribution (n=%d genes):\n", length(deg)))
  cat(sprintf("    Min:    %d\n", min(deg)))
  cat(sprintf("    Q1:     %d\n", as.integer(quantile(deg, 0.25))))
  cat(sprintf("    Median: %d\n", as.integer(median(deg))))
  cat(sprintf("    Q3:     %d\n", as.integer(quantile(deg, 0.75))))
  cat(sprintf("    Max:    %d\n", max(deg)))
  cat(sprintf("    Mean:   %.1f\n", mean(deg)))
}

# Community size summary
if (sum(!is.na(atlas$network_community_macro)) > 0) {
  macro_tbl <- table(atlas$network_community_macro)
  cat(sprintf("\n  Macro communities: %d communities\n", length(macro_tbl)))
  cat(sprintf("    Size range: %d - %d genes\n", min(macro_tbl), max(macro_tbl)))
  cat(sprintf("    Median size: %d genes\n", as.integer(median(macro_tbl))))
}

if (sum(!is.na(atlas$network_community_meso)) > 0) {
  meso_tbl <- table(atlas$network_community_meso)
  cat(sprintf("\n  Meso communities: %d communities\n", length(meso_tbl)))
  cat(sprintf("    Size range: %d - %d genes\n", min(meso_tbl), max(meso_tbl)))
  cat(sprintf("    Median size: %d genes\n", as.integer(median(meso_tbl))))
}

# Multiplicity distribution
if (sum(!is.na(atlas$network_edge_multiplicity_max)) > 0) {
  mult <- atlas$network_edge_multiplicity_max[!is.na(atlas$network_edge_multiplicity_max)]
  cat(sprintf("\n  Max edge multiplicity (n=%d genes):\n", length(mult)))
  cat(sprintf("    Distribution: %s\n",
              paste(names(table(mult)), table(mult), sep = "=", collapse = ", ")))
}

cat(sprintf("\nDone. Elapsed: %.1f seconds\n", as.numeric(difftime(Sys.time(),
    as.POSIXct(paste(Sys.Date(), "00:00:00")), units = "secs"))))
cat("End time:", format(Sys.time()), "\n")
