#!/usr/bin/env Rscript
# 34_combat_seq_sensitivity.R
# ---------------------------------------------------------------------------
# ComBat-seq sensitivity analysis for batch effects (Audit Issue 5)
#
# Variance partition showed dataset (batch) explains 22.7% of variance
# vs disease condition 0.91%. This script:
#   1. Applies ComBat-seq to remove batch effects
#   2. Re-runs dream on adjusted counts
#   3. Compares top DEGs with primary (unadjusted) results
#   4. Reports concordance metrics
#
# Output: RNA-seq/results/audit_sensitivity/combat_seq/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(sva)
  library(edgeR)
  library(variancePartition)
  library(ggplot2)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/combat_seq")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== ComBat-seq Sensitivity Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ================================================================
# 1. Load merged counts and metadata
# ================================================================
cat("--- Step 1: Loading data ---\n")

dge_file <- file.path(RDIR, "merged_dge.rds")
if (!file.exists(dge_file)) {
  cat("ERROR: merged_dge.rds not found. Run Script 03 first.\n")
  quit(save = "no", status = 1)
}

dge <- readRDS(dge_file)
cat("  Samples:", ncol(dge), "\n")
cat("  Genes:", nrow(dge), "\n")

meta <- as.data.table(dge$samples)
cat("  Datasets:", paste(unique(meta$dataset), collapse = ", "), "\n")

# Use group_binary (MASLD vs Control) to match primary dream analysis
if ("group_binary" %in% names(meta)) {
  condition_var <- "group_binary"
} else if ("condition" %in% names(meta)) {
  condition_var <- "condition"
} else {
  cat("  Available columns:", paste(names(meta), collapse = ", "), "\n")
  stop("Cannot find condition column")
}
cat("  Condition variable:", condition_var, "\n")
cat("  Condition levels:", paste(unique(meta[[condition_var]]), collapse = ", "), "\n")

# ================================================================
# 2. Apply ComBat-seq
# ================================================================
cat("\n--- Step 2: Applying ComBat-seq ---\n")

batch <- as.factor(meta$dataset)
group <- as.factor(meta[[condition_var]])

# ComBat-seq works on raw counts
counts_raw <- dge$counts

cat("  Running ComBat_seq (this may take a few minutes)...\n")
counts_adjusted <- tryCatch({
  ComBat_seq(counts = as.matrix(counts_raw),
             batch = batch,
             group = group)
}, error = function(e) {
  cat("  WARNING: ComBat_seq failed:", conditionMessage(e), "\n")
  cat("  Trying with covar_mod=NULL...\n")
  ComBat_seq(counts = as.matrix(counts_raw),
             batch = batch)
})

cat("  ComBat-seq complete. Adjusted counts dimensions:", dim(counts_adjusted), "\n")

# ================================================================
# 3. Re-run dream on adjusted counts
# ================================================================
cat("\n--- Step 3: Re-running dream on ComBat-seq adjusted counts ---\n")

# Create new DGEList with adjusted counts
dge_adj <- DGEList(counts = counts_adjusted, samples = meta)
dge_adj <- calcNormFactors(dge_adj, method = "TMM")

# Filter low-expression genes (same as primary analysis)
keep <- filterByExpr(dge_adj, group = meta[[condition_var]])
dge_adj <- dge_adj[keep, , keep.lib.sizes = FALSE]
cat("  Genes after filtering:", nrow(dge_adj), "\n")

# Voom
vobjDream <- voomWithDreamWeights(dge_adj,
                                   formula = as.formula(paste0("~ ", condition_var, " + (1|dataset)")),
                                   data = meta,
                                   BPPARAM = BiocParallel::MulticoreParam(4))

# Fit dream (still with random intercept for dataset — should be much smaller now)
cat("  Fitting dream model...\n")
fit <- dream(vobjDream,
             formula = as.formula(paste0("~ ", condition_var, " + (1|dataset)")),
             data = meta,
             BPPARAM = BiocParallel::MulticoreParam(4))

# Extract results for disease condition
coef_name <- grep("condition|disease", colnames(fit$coefficients), value = TRUE)
if (length(coef_name) == 0) {
  coef_name <- colnames(fit$coefficients)[2]  # Second coefficient is usually disease
}
cat("  Using coefficient:", coef_name[1], "\n")

tt <- topTable(fit, coef = coef_name[1], number = Inf, sort.by = "none")
tt_dt <- as.data.table(tt, keep.rownames = "gene")
cat("  Dream results:", nrow(tt_dt), "genes\n")
cat("  DEGs (padj<0.1):", sum(tt_dt$adj.P.Val < 0.1, na.rm = TRUE), "\n")

fwrite(tt_dt, file.path(OUTDIR, "dream_combat_seq_results.csv"))
cat("  Saved: dream_combat_seq_results.csv\n")

# ================================================================
# 4. Compare with primary (unadjusted) results
# ================================================================
cat("\n--- Step 4: Concordance with primary results ---\n")

primary_file <- file.path(RDIR, "dream_results.csv")
if (file.exists(primary_file)) {
  primary <- fread(primary_file)
  primary[, gene_clean := sub("\\..*", "", gene)]
  tt_dt[, gene_clean := sub("\\..*", "", gene)]

  # Primary uses "padj"; topTable output uses "adj.P.Val"
  padj_col_primary <- if ("padj" %in% names(primary)) "padj" else "adj.P.Val"
  padj_col_combat  <- if ("padj" %in% names(tt_dt)) "padj" else "adj.P.Val"

  merged <- merge(primary[, .(gene_clean, logFC_primary = logFC, padj_primary = get(padj_col_primary))],
                  tt_dt[, .(gene_clean, logFC_combat = logFC, padj_combat = get(padj_col_combat))],
                  by = "gene_clean")

  cat("  Genes in both analyses:", nrow(merged), "\n")

  # LFC correlation
  lfc_cor <- cor(merged$logFC_primary, merged$logFC_combat, use = "complete.obs")
  cat(sprintf("  logFC correlation (primary vs ComBat-seq): r = %.3f\n", lfc_cor))

  # Concordance of significant genes
  sig_primary <- merged[padj_primary < 0.1]$gene_clean
  sig_combat <- merged[padj_combat < 0.1]$gene_clean
  overlap <- length(intersect(sig_primary, sig_combat))
  jaccard <- overlap / length(union(sig_primary, sig_combat))

  cat(sprintf("  Significant genes — primary: %d, ComBat-seq: %d, overlap: %d\n",
              length(sig_primary), length(sig_combat), overlap))
  cat(sprintf("  Jaccard index: %.3f\n", jaccard))

  # Direction concordance for overlapping significant genes
  overlap_dt <- merged[gene_clean %in% intersect(sig_primary, sig_combat)]
  if (nrow(overlap_dt) > 0) {
    same_dir <- sum(sign(overlap_dt$logFC_primary) == sign(overlap_dt$logFC_combat))
    cat(sprintf("  Direction concordance in overlap: %d/%d (%.1f%%)\n",
                same_dir, nrow(overlap_dt), 100*same_dir/nrow(overlap_dt)))
  }

  # Save concordance metrics
  concordance_dt <- data.table(
    metric = c("n_genes_both", "lfc_pearson_r",
               "sig_primary", "sig_combat", "sig_overlap", "jaccard",
               "direction_concordance"),
    value = c(nrow(merged), round(lfc_cor, 4),
              length(sig_primary), length(sig_combat), overlap, round(jaccard, 4),
              if (nrow(overlap_dt) > 0) round(same_dir/nrow(overlap_dt), 4) else NA)
  )
  fwrite(concordance_dt, file.path(OUTDIR, "combat_seq_concordance.csv"))
  cat("  Saved: combat_seq_concordance.csv\n")

  # LFC scatter plot
  pdf(file.path(OUTDIR, "lfc_scatter_primary_vs_combat.pdf"), width = 6, height = 6)
  ggplot(merged, aes(x = logFC_primary, y = logFC_combat)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = sprintf("r = %.3f", lfc_cor), color = "red", size = 4) +
    labs(title = "logFC Concordance: Primary vs ComBat-seq",
         x = "logFC (primary dream)", y = "logFC (ComBat-seq + dream)") +
    theme_bw(base_size = 11)
  dev.off()
  cat("  Saved: lfc_scatter_primary_vs_combat.pdf\n")
} else {
  cat("  Primary dream results not found — cannot compare\n")
}

cat("\n=== ComBat-seq Sensitivity Analysis Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
