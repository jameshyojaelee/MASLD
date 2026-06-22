#!/usr/bin/env Rscript
# M04_mouse_consensus.R
# ---------------------------------------------------------------------------
# Mouse DEGs from the pooled limma-voom quality-weighted mega-analysis (M03;
# replaces dream 2026-06-16, mirroring the human limma_voom_qw__C2 canonical).
# Thresholds: padj < 0.05, |LFC| > 0.5  (cross-species match with human raw Tier-1).
# Input:  lvqw_pooled_results.csv
# Output: mouse_consensus_degs.csv
# NOTE: the emitted column names retain the legacy `dream_*` prefix purely for
#   downstream compatibility (the publication figure scripts/figures/figS06_cross_species.R
#   and retired Script 27 read those names). The VALUES are now limma-voom-qw, not dream.
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

cat("=== M04: Mouse pooled limma-voom-qw DEGs ===\n\n")

# --- Load results (limma-voom-qw pooled mega-analysis) ---
dream  <- fread(file.path(METADIR, "lvqw_pooled_results.csv"))

cat("Pooled limma-voom-qw results:", nrow(dream), "genes\n")

# --- Mark significance (legacy column names retained for downstream compatibility) ---
dream[, dream_logFC := logFC]
dream[, dream_padj := adj.P.Val]
dream[, dream_sig := !is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.5]  # raw Tier-1, matches human canonical |LFC|>0.5
dream[, direction := fcase(
  dream_logFC > 0, "Up",
  dream_logFC < 0, "Down",
  default = "Unknown"
)]

# Provide backwards compatibility wrapper column `tier` used by other scripts
# (treat significant as replacing Tier 1+2)
dream[, tier := ifelse(dream_sig, "Significant", "Not_significant")]

# --- Summary ---
cat("===== POOLED limma-voom-qw MEGA-ANALYSIS SUMMARY =====\n")
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
