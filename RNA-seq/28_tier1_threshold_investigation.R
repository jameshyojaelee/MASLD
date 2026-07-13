#!/usr/bin/env Rscript
# 28_tier1_threshold_investigation.R
# Diagnostic: why is tier1_high_confidence_degs.csv empty?

library(data.table)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
RDIR <- file.path(BASE, "analysis/integration/results/integration")

degs <- fread(file.path(RDIR, "consensus_degs.csv"))
cat("Total genes:", nrow(degs), "\n")
cat("Current tier distribution:\n")
print(table(degs$tier))

# Check bulk-only significance at various thresholds
for (padj_thresh in c(0.01, 0.05, 0.1)) {
  for (lfc_thresh in c(0.5, 0.8, 1.0)) {
    bulk_sig <- degs[bulk_padj < padj_thresh & abs(bulk_logFC) > lfc_thresh]
    meta_sig <- degs[meta_padj < padj_thresh & abs(meta_logFC) > lfc_thresh]
    both_sig <- degs[bulk_padj < padj_thresh & abs(bulk_logFC) > lfc_thresh &
                     meta_padj < padj_thresh & abs(meta_logFC) > lfc_thresh]
    # Check concordant direction
    concordant <- both_sig[sign(bulk_logFC) == sign(meta_logFC)]
    cat(sprintf("padj<%s |LFC|>%s: bulk=%d meta=%d both=%d concordant=%d\n",
                padj_thresh, lfc_thresh, nrow(bulk_sig), nrow(meta_sig),
                nrow(both_sig), nrow(concordant)))
  }
}

# Recommend: what threshold gives a reasonable Tier 1 count?
cat("\n--- Recommendation ---\n")
rec <- degs[bulk_padj < 0.05 & abs(bulk_logFC) > 0.5 &
            meta_padj < 0.05 & abs(meta_logFC) > 0.5 &
            sign(bulk_logFC) == sign(meta_logFC)]
cat("Recommended Tier 1 (padj<0.05, |LFC|>0.5, concordant):", nrow(rec), "genes\n")
if (nrow(rec) > 0) {
  cat("Top 20 by bulk padj:\n")
  print(head(rec[order(bulk_padj), .(gene, bulk_logFC, bulk_padj, meta_logFC, meta_padj)], 20))
}
