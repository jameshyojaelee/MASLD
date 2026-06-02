#!/usr/bin/env Rscript
# =============================================================================
# B6 Breakpoint sensitivity (V5 verification)
# Re-runs switch-like classification (Script 218) with breakpoints at F1, F2, F3.
# For each: linear model vs step model (fib >= k), per-gene AIC comparison.
# Compare AICs and step-like counts to verify F2 is the actual switch location.
#
# Pass: F2 has lowest mean AIC_step AND highest n_step_like among k∈{1,2,3}.
# Output: verification/controls/b6_breakpoint_sensitivity.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUT_CSV <- file.path(BASE, "docs/manuscript/verification/controls/b6_breakpoint_sensitivity.csv")
OUT_PER_GENE <- file.path(BASE, "docs/manuscript/verification/controls/b6_breakpoint_per_gene.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

cat("=== B6 breakpoint sensitivity ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ---- Load ----
dge <- readRDS(file.path(INT, "results/integration/merged_dge.rds"))
dream <- fread(file.path(INT, "results/integration/dream_results_ashr.csv"))
degs <- dream[padj < 0.05 & abs(logFC) > 0.3]
cat("DEGs:", nrow(degs), "\n")

meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"),
              select = c("sample_id", "fibrosis_stage", "dataset"))
meta <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta[, fibrosis_stage := as.integer(fibrosis_stage)]

shared_samples <- intersect(meta$sample_id, colnames(dge))
meta <- meta[sample_id %in% shared_samples]
shared_genes <- intersect(unique(degs$gene), rownames(dge))
cat("Shared samples:", length(shared_samples), " DEGs in DGE:", length(shared_genes), "\n")

dge_sub <- dge[shared_genes, shared_samples]
logcpm <- edgeR::cpm(dge_sub, log = TRUE, prior.count = 1)
stage_vec <- meta[match(colnames(logcpm), sample_id), fibrosis_stage]
fib_numeric <- as.numeric(stage_vec)

# ---- For each breakpoint, fit step model and compare to linear ----
breakpoints <- c(1L, 2L, 3L)
cat("\nFitting per-gene models for", nrow(logcpm), "genes × 3 breakpoints...\n")

# Linear model is fixed across breakpoints (fit once)
linear_aic <- numeric(nrow(logcpm))
for (i in seq_len(nrow(logcpm))) {
  fit_lin <- lm(logcpm[i, ] ~ fib_numeric)
  linear_aic[i] <- AIC(fit_lin)
}
cat("Linear AIC done. mean=", round(mean(linear_aic), 2), "\n")

# Step models at each breakpoint
per_gene <- data.table(gene = rownames(logcpm), aic_linear = linear_aic)

for (k in breakpoints) {
  cat(sprintf("  Breakpoint k=%d (fib >= %d) ...\n", k, k))
  fib_bin <- as.integer(stage_vec >= k)
  step_aic <- numeric(nrow(logcpm))
  for (i in seq_len(nrow(logcpm))) {
    fit_step <- lm(logcpm[i, ] ~ fib_bin)
    step_aic[i] <- AIC(fit_step)
    if (i %% 2000 == 0) cat(sprintf("    %d/%d\n", i, nrow(logcpm)))
  }
  per_gene[, paste0("aic_step_k", k) := step_aic]
  per_gene[, paste0("delta_aic_k", k) := aic_linear - step_aic]
}

fwrite(per_gene, OUT_PER_GENE)
cat("Per-gene saved:", OUT_PER_GENE, "\n")

# ---- Summarise ----
summary_rows <- lapply(breakpoints, function(k) {
  delta_col <- paste0("delta_aic_k", k)
  aic_col <- paste0("aic_step_k", k)
  n_step_gt0 <- sum(per_gene[[delta_col]] > 0, na.rm = TRUE)
  n_step_gt2 <- sum(per_gene[[delta_col]] > 2, na.rm = TRUE)
  n_step_gt4 <- sum(per_gene[[delta_col]] > 4, na.rm = TRUE)
  mean_aic_step <- mean(per_gene[[aic_col]], na.rm = TRUE)
  median_aic_step <- median(per_gene[[aic_col]], na.rm = TRUE)
  mean_delta <- mean(per_gene[[delta_col]], na.rm = TRUE)
  data.table(
    breakpoint = sprintf("F%d", k - 1L),  # i.e., k=2 → breakpoint "F1/F2" = between F1 and F2
    breakpoint_binary = sprintf("fib>=%d", k),
    n_genes = nrow(per_gene),
    n_step_like_delta_gt0 = n_step_gt0,
    pct_step_like_delta_gt0 = 100 * n_step_gt0 / nrow(per_gene),
    n_step_like_delta_gt2 = n_step_gt2,
    n_step_like_delta_gt4 = n_step_gt4,
    mean_aic_step = mean_aic_step,
    median_aic_step = median_aic_step,
    mean_delta_aic = mean_delta
  )
})
summary_dt <- rbindlist(summary_rows)
# Identify best breakpoint
summary_dt[, is_best_by_count := n_step_like_delta_gt0 == max(n_step_like_delta_gt0)]
summary_dt[, is_best_by_aic := mean_aic_step == min(mean_aic_step)]

cat("\n=== SUMMARY ===\n")
print(summary_dt)

fwrite(summary_dt, OUT_CSV)
cat("\nSaved:", OUT_CSV, "\n")
cat("End:", format(Sys.time()), "\n")
