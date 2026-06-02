#!/usr/bin/env Rscript
# aggregate_loo_cv.R
# ---------------------------------------------------------------------------
# Aggregate dream LOO-CV results: per-iteration metrics + per-gene stability.
# Run AFTER all 8 LOO-CV dream jobs complete.
# Output: loo_cv_summary.csv, loo_cv_per_gene.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR <- file.path(RDIR, "loo_cv")

# Cohorts mirror the mega-analysis set defined in config/human_datasets.yaml
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
COHORTS <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("LOO cohort set (k =", length(COHORTS), "):",
    paste(COHORTS, collapse = ", "), "\n")

# --- Load full model dream results ---
cat("Loading full dream results...\n")
full <- fread(file.path(RDIR, "dream_results.csv"))
full_degs <- full[padj < 0.1, gene]
cat(sprintf("Full model: %d genes tested, %d DEGs (padj<0.1)\n", nrow(full), length(full_degs)))

# --- Load all LOO results ---
cat("\nLoading LOO results...\n")
loo_list <- list()
missing <- character(0)

for (study in COHORTS) {
  f <- file.path(LOO_DIR, paste0("dream_loo_", study, ".csv"))
  if (!file.exists(f)) {
    cat(sprintf("  WARNING: Missing %s\n", f))
    missing <- c(missing, study)
    next
  }
  loo_list[[study]] <- fread(f)
  cat(sprintf("  %s: %d genes, %d DEGs\n", study, nrow(loo_list[[study]]),
              sum(loo_list[[study]]$padj < 0.1)))
}

if (length(missing) > 0) {
  cat(sprintf("\n*** %d LOO results missing: %s ***\n", length(missing), paste(missing, collapse = ", ")))
  cat("Run those LOO jobs first.\n")
  if (length(loo_list) == 0) stop("No LOO results found!")
}

# ===========================================================================
# Per-iteration metrics
# ===========================================================================
cat("\n── Computing per-iteration metrics ──\n")

iter_metrics <- rbindlist(lapply(names(loo_list), function(study) {
  loo <- loo_list[[study]]
  loo_degs <- loo[padj < 0.1, gene]

  # Merge on common genes
  merged <- merge(full[, .(gene, full_logFC = logFC, full_padj = padj)],
                  loo[, .(gene, loo_logFC = logFC, loo_padj = padj)],
                  by = "gene")

  # Spearman correlation of logFC
  rho <- cor(merged$full_logFC, merged$loo_logFC, method = "spearman")

  # Direction concordance (among genes significant in both)
  both_sig <- merged[full_padj < 0.1 & loo_padj < 0.1]
  dir_conc <- if (nrow(both_sig) > 0) {
    mean(sign(both_sig$full_logFC) == sign(both_sig$loo_logFC)) * 100
  } else NA_real_

  # Jaccard
  intersection <- length(intersect(full_degs, loo_degs))
  union_set <- length(union(full_degs, loo_degs))
  jaccard <- if (union_set > 0) intersection / union_set else 0

  # Recovery of full-model DEGs
  pct_recovered <- if (length(full_degs) > 0) {
    length(intersect(full_degs, loo_degs)) / length(full_degs) * 100
  } else 0

  data.table(
    held_out = study,
    n_degs = length(loo_degs),
    pct_full_recovered = round(pct_recovered, 1),
    spearman_rho = round(rho, 4),
    direction_concordance = round(dir_conc, 1),
    jaccard = round(jaccard, 4),
    n_genes_tested = nrow(loo)
  )
}))

cat("\n===== LOO-CV Summary =====\n")
print(iter_metrics)
cat(sprintf("\nMean DEGs: %.0f (full: %d)\n", mean(iter_metrics$n_degs), length(full_degs)))
cat(sprintf("Mean recovery: %.1f%%\n", mean(iter_metrics$pct_full_recovered)))
cat(sprintf("Mean Spearman rho: %.4f\n", mean(iter_metrics$spearman_rho)))
cat(sprintf("Mean direction concordance: %.1f%%\n", mean(iter_metrics$direction_concordance, na.rm = TRUE)))
cat(sprintf("Mean Jaccard: %.4f\n", mean(iter_metrics$jaccard)))

fwrite(iter_metrics, file.path(LOO_DIR, "loo_cv_summary.csv"))
cat("\nSaved: loo_cv_summary.csv\n")

# ===========================================================================
# Per-gene stability metrics (vectorized)
# ===========================================================================
cat("\n── Computing per-gene stability ──\n")

# Stack all LOO results into a single data.table for vectorized computation
all_loo <- rbindlist(loo_list, idcol = "held_out")

# Per-gene summary across LOO iterations
gene_loo_stats <- all_loo[, .(
  n_loo_sig    = sum(padj < 0.1),
  n_loo_tested = .N,
  min_logFC    = min(logFC),
  max_logFC    = max(logFC),
  mean_logFC   = mean(logFC),
  sd_logFC     = if (.N > 1) sd(logFC) else NA_real_
), by = gene]

# Merge with full-model significance
gene_stability <- merge(
  full[, .(gene, full_sig = padj < 0.1)],
  gene_loo_stats,
  by = "gene", all.x = TRUE
)

# Genes in full model but absent from all LOOs (rare; would have NA stats)
gene_stability[is.na(n_loo_tested), `:=`(n_loo_sig = 0L, n_loo_tested = 0L)]

# Robustness classification: relative to how many LOOs actually tested the gene
gene_stability[, robustness := fifelse(
  n_loo_tested == 0, "Untested",
  fifelse(n_loo_sig == n_loo_tested, "Robust",
  fifelse(n_loo_sig >= n_loo_tested - 2, "Stable",
          fifelse(n_loo_sig >= ceiling(n_loo_tested / 2), "Moderate", "Fragile")))
)]

# Summary
n_loo <- length(loo_list)
cat(sprintf("\nFull-model DEGs: %d\n", sum(gene_stability$full_sig)))
cat(sprintf("Robust (sig in all tested LOOs): %d\n",
            sum(gene_stability$full_sig & gene_stability$robustness == "Robust")))
cat(sprintf("Stable (sig in n-2+ tested LOOs): %d\n",
            sum(gene_stability$full_sig & gene_stability$robustness %in% c("Robust", "Stable"))))
cat(sprintf("Fragile (sig in <50%% tested LOOs): %d\n",
            sum(gene_stability$full_sig & gene_stability$robustness == "Fragile")))

fwrite(gene_stability, file.path(LOO_DIR, "loo_cv_per_gene.csv"))
cat("\nSaved: loo_cv_per_gene.csv\n")
cat("Done!\n")
