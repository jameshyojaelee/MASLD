#!/usr/bin/env Rscript
# M04_mouse_consensus.R
# ---------------------------------------------------------------------------
# Mouse DEGs from dream mega-analysis.
# Thresholds: padj < 0.05, |LFC| > 0.5  (2026-04-27 migration; cross-species match)
# Input:  dream_pooled_results.csv
# Output: mouse_consensus_degs.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
})

MOUSE  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse"
INT    <- file.path(MOUSE, "Unified_Integration")
RDIR   <- file.path(INT, "results")
METADIR <- file.path(RDIR, "meta_analysis")

cat("=== M04: Mouse Dream DEGs ===\n\n")

# --- Load results ---
dream  <- fread(file.path(METADIR, "dream_pooled_results.csv"))

cat("Dream results:", nrow(dream), "genes\n")

# --- Mark significance ---
dream[, dream_logFC := logFC]
dream[, dream_padj := adj.P.Val]
dream[, dream_sig := !is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.5]  # 2026-04-27: 1.0 -> 0.5 (match human canonical; mouse LOO-CV stability had no clean elbow)
dream[, direction := fcase(
  dream_logFC > 0, "Up",
  dream_logFC < 0, "Down",
  default = "Unknown"
)]

# Provide backwards compatibility wrapper column `tier` used by other scripts 
# (treat significant as replacing Tier 1+2)
dream[, tier := ifelse(dream_sig, "Significant", "Not_significant")]

# --- Summary ---
cat("===== DREAM MEGA-ANALYSIS SUMMARY =====\n")
tier_summary <- dream[, .N, by = tier][order(tier)]
print(tier_summary)

sig_genes <- dream[dream_sig == TRUE]
if (nrow(sig_genes) > 0) {
  cat("\nDirection (Significant):\n")
  print(table(sig_genes$tier, sig_genes$direction))
}

# --- Save ---
# Save with the historical name mouse_consensus_degs.csv for pipeline compatibility
fwrite(dream[order(-dream_sig, dream_padj)],
       file.path(RDIR, "mouse_consensus_degs.csv"))
cat("\nSaved:", file.path(RDIR, "mouse_consensus_degs.csv"), "\n")

cat("\nM04 complete.\n")
