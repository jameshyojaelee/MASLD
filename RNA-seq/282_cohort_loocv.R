#!/usr/bin/env Rscript
#SBATCH --job-name=net_282_loco
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_282_loco_%j.out
#SBATCH --error=logs/net_282_loco_%j.err
#
# Script 282: Leave-one-cohort-out (LOCO) replication for D-F2 edges.
# Implements falsification Criterion 1 (median replication fraction >= 0.70).
#
# Inputs:
#   - RNA-seq/results/network/edges_d_f2.csv          (from Script 281)
#   - merged_dge.rds                                  (per-sample log-CPM)
#   - RNA-seq/results/network/stage/stage_sample_manifest.csv
#
# Outputs:
#   - RNA-seq/results/network/edges_d_f2_loco.csv
#   - RNA-seq/results/network/stage/loco_summary.csv

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(matrixStats)
  library(parallel)
})

# ---- Paths ----------------------------------------------------------------
PROJ       <- Sys.getenv("MASLD_PROJECT_ROOT",
                         "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
NET_DIR    <- file.path(PROJ, "RNA-seq/results/network")
STAGE_DIR  <- file.path(NET_DIR, "stage")
DGE_PATH   <- file.path(PROJ, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
EDGE_PATH  <- file.path(NET_DIR, "edges_d_f2.csv")
MANIFEST   <- file.path(STAGE_DIR, "stage_sample_manifest.csv")

OUT_EDGES   <- file.path(NET_DIR, "edges_d_f2_loco.csv")
OUT_SUMMARY <- file.path(STAGE_DIR, "loco_summary.csv")

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))

# ---- Validate inputs ------------------------------------------------------
stopifnot(file.exists(EDGE_PATH))
stopifnot(file.exists(DGE_PATH))
stopifnot(file.exists(MANIFEST))

# ---- Load -----------------------------------------------------------------
cat("[282] Loading D-F2 edges, DGE, manifest...\n")
edges    <- fread(EDGE_PATH)
dge      <- readRDS(DGE_PATH)

# Build per-sample manifest directly from unified_metadata.csv (cohort + fibrosis_stage)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
UNIFIED_META <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
unified_meta <- fread(UNIFIED_META)

norm_stage <- function(x) {
  x <- as.character(x)
  x[x %in% c("0", "1", "F0", "F1")] <- "F01"
  x[x %in% c("2", "F2")]             <- "F2"
  x[x %in% c("3", "4", "F3", "F4")]  <- "F34"
  x[!(x %in% c("F01", "F2", "F34"))] <- NA_character_
  x
}

manifest <- data.table(sample_id = colnames(dge))
manifest <- merge(manifest, unified_meta[, .(sample_id, dataset, fibrosis_stage)],
                  by = "sample_id", all.x = TRUE, sort = FALSE)
setnames(manifest, "dataset", "cohort")
manifest[, stage_bin := norm_stage(fibrosis_stage)]
manifest <- manifest[!is.na(cohort) & !is.na(stage_bin)]
cohorts <- sort(unique(manifest$cohort))
cat(sprintf("[282] %d cohorts: %s\n", length(cohorts), paste(cohorts, collapse = ", ")))

# ---- Compute logCPM -------------------------------------------------------
cat("[282] Computing logCPM...\n")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# ENSG→symbol mapping from node set (same pattern as Script 280)
NODE_PATH <- file.path(BASE, "RNA-seq/results/network/network_nodes.csv")
nodes <- fread(NODE_PATH)
strip_version <- function(ids) sub("\\.[0-9]+$", "", ids)
sym_map <- setNames(nodes$human_symbol, strip_version(nodes$ensembl_id))

rn_stripped <- strip_version(rownames(logcpm))
has_sym <- !is.na(sym_map[rn_stripped]) & nzchar(sym_map[rn_stripped])
new_rn <- ifelse(has_sym, sym_map[rn_stripped], rownames(logcpm))
# collapse duplicates (keep row with highest mean)
if (any(duplicated(new_rn[has_sym]))) {
  ord <- order(-rowMeans(logcpm, na.rm = TRUE))
  logcpm <- logcpm[ord, , drop = FALSE]
  new_rn <- new_rn[ord]
  keep_row <- !duplicated(new_rn)
  logcpm <- logcpm[keep_row, , drop = FALSE]
  new_rn <- new_rn[keep_row]
}
rownames(logcpm) <- new_rn

needed_genes <- unique(c(edges$gene_a, edges$gene_b))
gene_rn <- rownames(logcpm)
keep_genes <- intersect(needed_genes, gene_rn)
cat(sprintf("[282] %d / %d D-F2 edge genes found in DGE\n",
            length(keep_genes), length(needed_genes)))
logcpm <- logcpm[keep_genes, , drop = FALSE]

# Subset edges to those with both endpoints available
edges[, both_available := (gene_a %in% keep_genes) & (gene_b %in% keep_genes)]
edges_usable <- edges[both_available == TRUE]
cat(sprintf("[282] %d / %d edges usable\n", nrow(edges_usable), nrow(edges)))

# ---- Helper: compute per-stage correlation for a sample subset ------------
compute_delta_max <- function(sample_ids_by_stage, pair_idx_a, pair_idx_b) {
  # Returns data.table with cor_F01, cor_F2, cor_F34, delta_max
  stages <- c("F01", "F2", "F34")
  cors <- matrix(NA_real_, nrow = length(pair_idx_a), ncol = length(stages),
                 dimnames = list(NULL, stages))
  for (st in stages) {
    sids <- sample_ids_by_stage[[st]]
    sids <- intersect(sids, colnames(logcpm))
    if (length(sids) < 5) next
    M <- logcpm[, sids, drop = FALSE]
    # Compute pairwise correlations only for the needed rows; use per-pair loop
    va <- M[pair_idx_a, , drop = FALSE]
    vb <- M[pair_idx_b, , drop = FALSE]
    # Row-wise Pearson correlation between va[i,] and vb[i,]
    va_c <- va - rowMeans(va)
    vb_c <- vb - rowMeans(vb)
    num <- rowSums(va_c * vb_c)
    den <- sqrt(rowSums(va_c^2) * rowSums(vb_c^2))
    cors[, st] <- ifelse(den > 0, num / den, NA_real_)
  }
  # delta_max = max deviation of F2 from the mean of F01 and F34 (or max pairwise)
  d1 <- cors[, "F2"] - cors[, "F01"]
  d2 <- cors[, "F2"] - cors[, "F34"]
  # Choose the deviation with larger absolute value (same sign signal)
  pick <- ifelse(abs(d1) >= abs(d2) | is.na(d2), d1, d2)
  data.table(cor_F01 = cors[, "F01"],
             cor_F2  = cors[, "F2"],
             cor_F34 = cors[, "F34"],
             delta_max = pick)
}

# ---- Stage-grouped sample lists -------------------------------------------
stage_samples_all <- split(manifest$sample_id, manifest$stage_bin)

# Row indices for edge endpoints (into logcpm)
idx_a <- match(edges_usable$gene_a, keep_genes)
idx_b <- match(edges_usable$gene_b, keep_genes)

# ---- Full-cohort baseline delta_max ---------------------------------------
cat("[282] Computing full baseline delta_max...\n")
full_res <- compute_delta_max(stage_samples_all, idx_a, idx_b)
edges_usable[, delta_max_full := full_res$delta_max]

# If the upstream D-F2 table already has delta_max, we trust it; otherwise use
# our recomputation. Use the column that exists.
if (!"delta_max" %in% colnames(edges_usable)) {
  edges_usable[, delta_max := delta_max_full]
}

# ---- LOCO iterations ------------------------------------------------------
cat("[282] Running LOCO across cohorts...\n")
loco_mat <- matrix(0L, nrow = nrow(edges_usable), ncol = length(cohorts),
                   dimnames = list(NULL, cohorts))
per_cohort_rate <- numeric(length(cohorts))
names(per_cohort_rate) <- cohorts

for (k in seq_along(cohorts)) {
  C <- cohorts[k]
  keep_samples <- manifest[cohort != C]
  stage_samples_loco <- split(keep_samples$sample_id, keep_samples$stage_bin)
  res <- compute_delta_max(stage_samples_loco, idx_a, idx_b)
  d_loco <- res$delta_max
  d_full <- edges_usable$delta_max
  ok <- !is.na(d_loco) & !is.na(d_full) &
        (sign(d_loco) == sign(d_full)) &
        (abs(d_loco) >= 0.5 * abs(d_full))
  loco_mat[, k] <- as.integer(ok)
  per_cohort_rate[k] <- mean(ok, na.rm = TRUE)
  cat(sprintf("  [%2d/%d] drop %-25s  replicate_rate = %.3f\n",
              k, length(cohorts), C, per_cohort_rate[k]))
}

# ---- Aggregate ------------------------------------------------------------
replicates_count    <- rowSums(loco_mat)
replicates_fraction <- replicates_count / length(cohorts)

edges_usable[, replicates_loco_count    := replicates_count]
edges_usable[, replicates_loco_fraction := replicates_fraction]
edges_usable[, criterion1_status := ifelse(replicates_loco_fraction >= 0.70,
                                           "PASS", "FAIL")]

# Merge back onto full edges table (unusable edges get NA)
edges_merge <- merge(edges, edges_usable[, .(gene_a, gene_b,
                                             replicates_loco_count,
                                             replicates_loco_fraction,
                                             criterion1_status)],
                     by = c("gene_a", "gene_b"), all.x = TRUE, sort = FALSE)

fwrite(edges_merge, OUT_EDGES)
cat(sprintf("[282] Wrote %s (%d rows)\n", OUT_EDGES, nrow(edges_merge)))

# ---- Summary --------------------------------------------------------------
med_frac  <- median(replicates_fraction, na.rm = TRUE)
mean_frac <- mean(replicates_fraction, na.rm = TRUE)
pass_rate <- mean(replicates_fraction >= 0.70, na.rm = TRUE)

summary_dt <- rbind(
  data.table(metric = "per_cohort_replication_rate",
             cohort = names(per_cohort_rate),
             value  = per_cohort_rate),
  data.table(metric = "median_replication_fraction",  cohort = "ALL", value = med_frac),
  data.table(metric = "mean_replication_fraction",    cohort = "ALL", value = mean_frac),
  data.table(metric = "edges_passing_0.70_fraction",  cohort = "ALL", value = pass_rate),
  data.table(metric = "n_edges_tested",               cohort = "ALL",
             value = nrow(edges_usable))
)
dir.create(STAGE_DIR, showWarnings = FALSE, recursive = TRUE)
fwrite(summary_dt, OUT_SUMMARY)

cat("\n==================== CRITERION 1 ====================\n")
cat(sprintf("  median(replication_fraction) = %.3f  (threshold: 0.70)\n", med_frac))
cat(sprintf("  mean(replication_fraction)   = %.3f\n", mean_frac))
cat(sprintf("  edges passing >= 0.70        = %.1f%%\n", 100 * pass_rate))
cat(sprintf("  CRITERION 1: %s\n", if (med_frac >= 0.70) "PASS" else "FAIL"))
cat("=====================================================\n")
