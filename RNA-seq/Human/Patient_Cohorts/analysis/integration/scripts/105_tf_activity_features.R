#!/usr/bin/env Rscript
# 105_tf_activity_features.R
# ---------------------------------------------------------------------------
# Compute per-sample TF activity scores using decoupleR ULM on logCPM with
# DoRothEA regulons. The existing TF activity results in
# RNA-seq/results/multi_evidence/functional_activity/ are gene-level summaries
# (from dream t-statistics, 1 condition), not per-sample. This script computes
# the full samples x TFs activity matrix needed for the foundation model.
#
# Input:  results/integration/merged_dge.rds (34,453 genes x 1,444 samples)
# Output: results/staging_classifier/tf_activity_features.rds (samples x TFs matrix)
#         results/staging_classifier/tf_activity_names.csv
#
# Usage: Rscript 105_tf_activity_features.R
# SLURM: cpu, 8 CPUs, 64G, 48h
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(decoupleR)
})

set.seed(42)

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
ME     <- file.path(BASE, "RNA-seq/results/multi_evidence")
FA     <- file.path(ME, "functional_activity")
OUTDIR <- file.path(INT, "results/staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 105: Per-Sample TF Activity Features (decoupleR ULM) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# STEP 1: Load expression data and compute TMM logCPM
# ============================================================
cat("--- Step 1: Load DGE and compute logCPM ---\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("DGE loaded:", nrow(dge), "genes x", ncol(dge), "samples\n")

dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat("logCPM matrix:", nrow(logcpm), "x", ncol(logcpm), "\n\n")

# ============================================================
# STEP 2: Map Ensembl IDs to gene symbols
# ============================================================
cat("--- Step 2: Map Ensembl IDs to gene symbols ---\n")

atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
atlas <- atlas[!is.na(human_symbol) & human_symbol != ""]
atlas[, ensembl_clean := sub("\\..*", "", ensembl_id)]
atlas <- atlas[!duplicated(ensembl_clean)]

ensembl_clean <- sub("\\..*", "", rownames(logcpm))
idx <- match(ensembl_clean, atlas$ensembl_clean)
mapped <- !is.na(idx)
cat("Genes mapped to symbols:", sum(mapped), "of", nrow(logcpm),
    sprintf("(%.1f%%)\n", 100 * sum(mapped) / nrow(logcpm)))

logcpm_sym <- logcpm[mapped, ]
symbols <- atlas$human_symbol[idx[mapped]]

# Handle duplicate symbols: keep highest mean expression
mean_expr <- rowMeans(logcpm_sym)
dup_syms <- symbols[duplicated(symbols)]
if (length(dup_syms) > 0) {
  cat("Resolving", length(unique(dup_syms)), "duplicate symbols\n")
  keep <- rep(TRUE, length(symbols))
  for (s in unique(dup_syms)) {
    which_dup <- which(symbols == s)
    best <- which_dup[which.max(mean_expr[which_dup])]
    keep[setdiff(which_dup, best)] <- FALSE
  }
  logcpm_sym <- logcpm_sym[keep, ]
  symbols <- symbols[keep]
}

rownames(logcpm_sym) <- symbols
cat("Final symbol-mapped matrix:", nrow(logcpm_sym), "genes x",
    ncol(logcpm_sym), "samples\n\n")

# ============================================================
# STEP 3: Load DoRothEA TF-target network
# ============================================================
cat("--- Step 3: Load DoRothEA regulons (levels A, B, C) ---\n")
tf_net <- as.data.table(get_dorothea(organism = "human", levels = c("A", "B", "C")))
cat("DoRothEA network:", nrow(tf_net), "interactions,",
    length(unique(tf_net$source)), "TFs\n")

# Check overlap with our gene universe
net_targets <- unique(tf_net$target)
n_overlap <- sum(net_targets %in% rownames(logcpm_sym))
cat("Network targets in our data:", n_overlap, "of", length(net_targets),
    sprintf("(%.1f%%)\n\n", 100 * n_overlap / length(net_targets)))

# ============================================================
# STEP 4: Run decoupleR ULM per sample
# ============================================================
cat("--- Step 4: Run ULM on", ncol(logcpm_sym), "samples ---\n")
cat("This may take 10-30 minutes...\n")
t0 <- Sys.time()

# decoupleR run_ulm: mat = features (rows) x samples (cols), net = TF-target
# Returns long-format table with source (TF), condition (sample), score, p_value
ulm_results <- as.data.table(run_ulm(
  mat = logcpm_sym,
  net = tf_net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
))

t1 <- Sys.time()
cat("ULM completed in", round(difftime(t1, t0, units = "mins"), 1), "minutes\n")
cat("Raw results:", nrow(ulm_results), "rows\n")
cat("  TFs:", length(unique(ulm_results$source)), "\n")
cat("  Samples:", length(unique(ulm_results$condition)), "\n\n")

# ============================================================
# STEP 5: Reshape to samples x TFs matrix
# ============================================================
cat("--- Step 5: Reshape to activity matrix ---\n")

# Pivot: rows = samples, columns = TFs, values = ULM score
tf_wide <- dcast(ulm_results, condition ~ source, value.var = "score")
sample_ids <- tf_wide$condition
tf_wide[, condition := NULL]

tf_matrix <- as.matrix(tf_wide)
rownames(tf_matrix) <- sample_ids

# Reorder rows to match original sample order
sample_order <- colnames(logcpm)
tf_matrix <- tf_matrix[match(sample_order, rownames(tf_matrix)), ]
stopifnot(all(rownames(tf_matrix) == sample_order))

cat("TF activity matrix:", nrow(tf_matrix), "samples x", ncol(tf_matrix), "TFs\n")
cat("Score range:", round(range(tf_matrix, na.rm = TRUE), 3), "\n\n")

# ============================================================
# STEP 6: Cross-reference with existing gene-level TF results
# ============================================================
cat("--- Step 6: Cross-reference with existing results ---\n")

existing_tf <- file.path(FA, "tf_activity_scores.csv")
if (file.exists(existing_tf)) {
  existing <- fread(existing_tf)
  sig_tfs <- existing[padj < 0.05, tf]
  cat("Existing gene-level analysis: ", nrow(existing), " TFs tested, ",
      length(sig_tfs), " significant (padj < 0.05)\n", sep = "")

  # Check overlap
  our_tfs <- colnames(tf_matrix)
  cat("Our per-sample TFs:", length(our_tfs), "\n")
  cat("Overlap with significant TFs:", sum(sig_tfs %in% our_tfs),
      "of", length(sig_tfs), "\n")

  # Report top TFs from existing analysis that are in our matrix
  top_existing <- existing[1:10]
  cat("\nTop 10 TFs from existing analysis (gene-level dream t-stat):\n")
  for (i in seq_len(nrow(top_existing))) {
    tf_name <- top_existing$tf[i]
    in_ours <- tf_name %in% our_tfs
    cat(sprintf("  %2d. %s (score=%.2f, padj=%.1e) %s\n",
                i, tf_name, top_existing$score[i], top_existing$padj[i],
                if (in_ours) "[in matrix]" else "[NOT in matrix]"))
  }
} else {
  cat("No existing TF activity file found at:", existing_tf, "\n")
  sig_tfs <- character(0)
}

# ============================================================
# STEP 7: Filter to significant TFs (optional subset)
# ============================================================
cat("\n--- Step 7: Create filtered TF subset ---\n")

# Keep all TFs in the full matrix, but also create a filtered version
# with only the TFs significant in the gene-level analysis
if (length(sig_tfs) > 0) {
  sig_in_matrix <- intersect(sig_tfs, colnames(tf_matrix))
  tf_matrix_sig <- tf_matrix[, sig_in_matrix, drop = FALSE]
  cat("Significant-TF subset:", nrow(tf_matrix_sig), "samples x",
      ncol(tf_matrix_sig), "TFs\n")
} else {
  # If no existing analysis, keep top 70 TFs by variance across samples
  tf_vars <- apply(tf_matrix, 2, var)
  top70 <- names(sort(tf_vars, decreasing = TRUE))[1:min(70, ncol(tf_matrix))]
  tf_matrix_sig <- tf_matrix[, top70, drop = FALSE]
  cat("Top-70-by-variance subset:", nrow(tf_matrix_sig), "samples x",
      ncol(tf_matrix_sig), "TFs\n")
}

# ============================================================
# STEP 8: Save outputs
# ============================================================
cat("\n--- Step 8: Save outputs ---\n")

# Full TF activity matrix
saveRDS(tf_matrix, file.path(OUTDIR, "tf_activity_features.rds"))
cat("Saved:", file.path(OUTDIR, "tf_activity_features.rds"), "\n")

# Filtered TF activity matrix
saveRDS(tf_matrix_sig, file.path(OUTDIR, "tf_activity_features_significant.rds"))
cat("Saved:", file.path(OUTDIR, "tf_activity_features_significant.rds"), "\n")

# TF feature names
tf_names <- data.table(
  tf = colnames(tf_matrix),
  mean_activity = colMeans(tf_matrix),
  sd_activity = apply(tf_matrix, 2, sd),
  is_significant = colnames(tf_matrix) %in% sig_tfs
)
setorder(tf_names, -sd_activity)
fwrite(tf_names, file.path(OUTDIR, "tf_activity_names.csv"))
cat("Saved:", file.path(OUTDIR, "tf_activity_names.csv"), "\n")

cat("\n=== Summary ===\n")
cat("Full TF matrix:", nrow(tf_matrix), "samples x", ncol(tf_matrix), "TFs\n")
cat("Filtered TF matrix:", nrow(tf_matrix_sig), "samples x", ncol(tf_matrix_sig), "TFs\n")
cat("TFs with highest variance across samples:\n")
print(tf_names[1:10, .(tf, mean_activity = round(mean_activity, 3),
                         sd_activity = round(sd_activity, 3), is_significant)])
cat("\nFinished:", as.character(Sys.time()), "\n")
