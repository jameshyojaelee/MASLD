#!/usr/bin/env Rscript
# 07_consensus_degs.R
# ---------------------------------------------------------------------------
# Identify primary disease-associated DEGs from dream mega-analysis.
# Uses standard padj + logFC thresholds for DEG definition:
#   padj < 0.05 AND |logFC| > LFC_THRESH
# ashr columns (lfsr, shrunk_logFC) preserved for supplementary reference.
#
# Output: results/integration/consensus_degs.csv
#         results/integration/dream_significant_degs.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Load results (prefer ashr-shrunk) ---
ashr_file <- file.path(RDIR, "dream_results_ashr.csv")
raw_file  <- file.path(RDIR, "dream_results.csv")

if (file.exists(ashr_file)) {
  dream <- fread(ashr_file)
  has_ashr <- "shrunk_logFC" %in% names(dream) & "lfsr" %in% names(dream)
  cat("Loaded ashr-shrunk dream results:", nrow(dream), "genes\n")
} else {
  dream <- fread(raw_file)
  has_ashr <- FALSE
  cat("WARNING: ashr output not found — using raw dream results\n")
  cat("Run 05b_ashr_shrinkage.R first for principled thresholding\n")
}

# --- Define significance thresholds ---
if (has_ashr) {
  PADJ_THRESH <- 0.05
  LFC_THRESH <- 0.5  # Migrated 0.3 -> 0.5 (LOO-CV stability; ~1.41x fold change)

  dream[, dream_sig := padj < PADJ_THRESH & abs(logFC) > LFC_THRESH]
  dream[, dream_dir := sign(logFC)]

  cat("\n===== DREAM MEGA-ANALYSIS DE RESULTS (padj + logFC) =====\n")
  cat("Threshold: padj <", PADJ_THRESH, "AND |logFC| >", LFC_THRESH, "\n")
} else {
  # Fallback: original thresholds
  PADJ_THRESH <- 0.1
  LFC_THRESH <- 0.5

  dream[, dream_sig := padj < PADJ_THRESH & abs(logFC) > LFC_THRESH]
  dream[, dream_dir := sign(logFC)]

  cat("\n===== DREAM MEGA-ANALYSIS DE RESULTS (raw, no ashr) =====\n")
  cat("Threshold: padj <", PADJ_THRESH, "AND |logFC| >", LFC_THRESH, "\n")
}

cat("Total genes:", nrow(dream), "\n\n")

# --- Summary ---
sig_genes <- dream[dream_sig == TRUE]
cat("Significant DEGs:", nrow(sig_genes), "\n")
cat("  Up:", sum(sig_genes$dream_dir == 1, na.rm = TRUE), "\n")
cat("  Down:", sum(sig_genes$dream_dir == -1, na.rm = TRUE), "\n")

# --- Rename columns for downstream compatibility ---
# dream_logFC = raw logFC (primary effect size)
# dream_shrunk_logFC = ashr shrunk_logFC (supplementary reference)
# dream_padj = padj
# dream_lfsr = lfsr (supplementary reference)

if (has_ashr) {
  setnames(dream, "logFC", "dream_logFC")
  setnames(dream, "shrunk_logFC", "dream_shrunk_logFC")
  if ("padj" %in% names(dream)) setnames(dream, "padj", "dream_padj")
  if ("lfsr" %in% names(dream)) setnames(dream, "lfsr", "dream_lfsr")
} else {
  setnames(dream, c("logFC", "padj"), c("dream_logFC", "dream_padj"))
}

# --- Save ---
dream <- dream[order(dream$dream_padj)]
fwrite(dream, file.path(RDIR, "consensus_degs.csv"))
cat("\nSaved: consensus_degs.csv\n")

fwrite(dream[dream_sig == TRUE], file.path(RDIR, "dream_significant_degs.csv"))
cat("Saved: dream_significant_degs.csv (", nrow(sig_genes), "genes)\n")

cat("\n===== COMPLETE =====\n")
