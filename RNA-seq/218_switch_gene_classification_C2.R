#!/usr/bin/env Rscript
# =============================================================================
# 218_switch_gene_classification_C2.R
# Item 3 of the C2 validation refresh: re-run the switch-like-vs-gradual
# AIC-at-F2-boundary DEG classification on the C2 canonical DEG set.
#
# OLD (dream-era): 1,992 switch-like of 5,484 dream DEGs (padj<.05 & |logFC|>.3).
# NEW (C2): denominator = C2 canonical Tier-1 raw (padj<.05 & |logFC|>.5) = 1,853.
#
# The AIC classification itself is DE-method-invariant: it fits per-gene
# expression~stage (linear) vs expression~I(stage>=2) (step) on logCPM from
# merged_dge.rds and the fibrosis staging metadata — IDENTICAL machinery to
# script 218. Only the DEG LIST being classified changes (dream -> C2 canonical).
#
# Output:
#   RNA-seq/results/stratified_causal/switch_gene_classification_C2.csv
# (No figure — numbers-only refresh; figure regen lives in script 218.)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUT_CSV <- file.path(BASE, "RNA-seq/results/stratified_causal/switch_gene_classification_C2.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

cat("=== 218 (C2): Switch-like gene classification on C2 canonical DEGs ===\n")

# -- 1. Load DGE + C2 canonical DEG list (Tier-1 raw) ------------------------
dge <- readRDS(file.path(INT, "results/integration/merged_dge.rds"))
cat("DGE:", nrow(dge), "genes x", ncol(dge), "samples\n")

can <- fread(file.path(INT, "results/integration/canonical_deg_results.csv"))
# C2 Tier-1 raw threshold (matches 05h producer + project canonical):
degs <- can[!is.na(padj) & padj < 0.05 & abs(logFC) > 0.5]
cat("C2 canonical Tier-1 DEGs (padj<.05 & |logFC|>.5):", nrow(degs), "\n")
deg_genes <- unique(degs$gene)   # VERSIONED ENSG, matching merged_dge rownames

# -- 2. Fibrosis staging metadata -------------------------------------------
meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"),
              select = c("sample_id", "fibrosis_stage", "dataset"))
meta <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta[, fibrosis_stage := as.integer(fibrosis_stage)]
cat("Samples with fibrosis stage:", nrow(meta), "\n")
print(table(meta$fibrosis_stage))

# -- 3. logCPM matrix (DEGs x staged samples) -------------------------------
shared_samples <- intersect(meta$sample_id, colnames(dge))
meta <- meta[sample_id %in% shared_samples]
shared_genes <- intersect(deg_genes, rownames(dge))
cat("Shared samples:", length(shared_samples), " | DEGs in DGE:", length(shared_genes), "\n")

dge_sub <- dge[shared_genes, shared_samples]
logcpm  <- edgeR::cpm(dge_sub, log = TRUE, prior.count = 1)

# -- 4. Per-stage means + per-gene linear vs step AIC -----------------------
stage_vec <- meta[match(colnames(logcpm), sample_id), fibrosis_stage]
stage_idx <- lapply(0:4, function(s) which(stage_vec == s))
stage_means <- sapply(0:4, function(s) rowMeans(logcpm[, stage_idx[[s + 1]], drop = FALSE]))
colnames(stage_means) <- paste0("F", 0:4)

fib_numeric <- as.numeric(stage_vec)
fib_binary  <- as.integer(stage_vec >= 2)   # F0-F1 vs F2-F4 (the F2 boundary)
gene_names  <- rownames(logcpm)

results_list <- vector("list", nrow(logcpm))
for (i in seq_len(nrow(logcpm))) {
  y <- logcpm[i, ]
  fit_lin  <- lm(y ~ fib_numeric)
  fit_step <- lm(y ~ fib_binary)
  delta_aic <- AIC(fit_lin) - AIC(fit_step)   # positive = step (switch) wins
  m <- stage_means[i, ]
  mean_early <- mean(m[1:2]); mean_late <- mean(m[3:5])
  epsilon <- 0.01
  switch_ratio <- abs(mean_late - mean_early) /
    (abs(m[2] - m[1]) + abs(m[5] - m[4]) + epsilon)
  direction <- ifelse(mean_late > mean_early, "up", "down")
  results_list[[i]] <- data.table(
    gene = gene_names[i],
    mean_F0 = m[1], mean_F1 = m[2], mean_F2 = m[3], mean_F3 = m[4], mean_F4 = m[5],
    mean_early = mean_early, mean_late = mean_late,
    delta_aic = delta_aic, r2_linear = summary(fit_lin)$r.squared,
    r2_step = summary(fit_step)$r.squared, switch_ratio = switch_ratio,
    direction = direction)
  if (i %% 500 == 0) cat("  ", i, "/", nrow(logcpm), "\n")
}
res <- rbindlist(results_list)

# -- 5. Classify (delta_aic > 0 => switch-like) -----------------------------
res[, classification := fifelse(
  delta_aic > 0 & direction == "up", "switch_up",
  fifelse(delta_aic > 0 & direction == "down", "switch_down",
  fifelse(delta_aic <= 0 & direction == "up", "gradual_up", "gradual_down")))]
res <- merge(res, degs[, .(gene, logFC, padj, t)], by = "gene", all.x = TRUE)

fwrite(res[order(-switch_ratio)], OUT_CSV)

# -- 6. Summary -------------------------------------------------------------
n_switch  <- sum(grepl("switch",  res$classification))
n_gradual <- sum(grepl("gradual", res$classification))
cat("\n=== SUMMARY (C2) ===\n")
cat(sprintf("Denominator (DEGs classified): %d\n", nrow(res)))
cat(sprintf("Switch-like: %d (%.1f%%)\n", n_switch, 100 * n_switch / nrow(res)))
cat(sprintf("Gradual:     %d (%.1f%%)\n", n_gradual, 100 * n_gradual / nrow(res)))
cat("\nClassification table:\n"); print(table(res$classification))
cat(sprintf("\nMedian switch ratio: %.3f\n", median(res$switch_ratio)))
cat(sprintf("Mean delta AIC (step - linear): %.2f\n", mean(res$delta_aic)))
cat("\nSaved:", OUT_CSV, "\n")
