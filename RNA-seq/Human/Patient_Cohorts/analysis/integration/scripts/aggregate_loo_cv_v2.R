#!/usr/bin/env Rscript
# aggregate_loo_cv_v2.R
# ---------------------------------------------------------------------------
# Aggregate dream LOO-CV v2 results.
# Run AFTER all 5 per-fold LOO-CV v2 jobs complete.
#
# Outputs:
#   loo_cv_v2/loo_cv_summary.csv           — per-fold metrics + means
#   loo_cv_v2/loo_cv_per_gene_stability.csv — per-gene: in how many folds is
#                                             the gene significant?
#
# Also computes v1-style "stability" metrics (% of full-model DEGs recovered)
# alongside the new v2 "replication" metrics, so both can be reported.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR <- file.path(RDIR, "loo_cv_v2")

# Cohorts from yaml
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
COHORTS <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("LOO cohort set (k =", length(COHORTS), "):",
    paste(COHORTS, collapse = ", "), "\n")

# ===========================================================================
# 1. Collect per-fold metrics
# ===========================================================================
cat("\n── Collecting per-fold metrics ──\n")

metrics_list <- list()
missing <- character(0)

for (study in COHORTS) {
  f <- file.path(LOO_DIR, paste0("loo_v2_metrics_", study, ".csv"))
  if (!file.exists(f)) {
    cat(sprintf("  WARNING: Missing %s\n", f))
    missing <- c(missing, study)
    next
  }
  metrics_list[[study]] <- fread(f)
  cat(sprintf("  %s: loaded\n", study))
}

if (length(missing) > 0) {
  cat(sprintf("\n*** %d fold results missing: %s ***\n",
              length(missing), paste(missing, collapse = ", ")))
  if (length(metrics_list) == 0) stop("No LOO-CV v2 results found!")
}

summary_dt <- rbindlist(metrics_list, fill = TRUE)

# ===========================================================================
# 2. Add v1-style stability metrics (training vs full model)
# ===========================================================================
cat("\n── Adding training-vs-full stability metrics ──\n")

full_file <- file.path(RDIR, "dream_results.csv")
if (file.exists(full_file)) {
  full <- fread(full_file)
  full_degs_01  <- full[padj < 0.1, gene]
  full_degs_005 <- full[padj < 0.05, gene]
  full_degs_lfc <- full[padj < 0.05 & abs(logFC) >= 0.5, gene]
  cat(sprintf("Full model: %d DEGs (padj<0.1), %d (padj<0.05), %d (padj<0.05+|lfc|>0.5)\n",
              length(full_degs_01), length(full_degs_005), length(full_degs_lfc)))

  # For each fold, compute stability = % of full-model DEGs recovered
  stability_cols <- rbindlist(lapply(names(metrics_list), function(study) {
    loo_file <- file.path(LOO_DIR, paste0("dream_loo_v2_", study, ".csv"))
    if (!file.exists(loo_file)) return(data.table(held_out_cohort = study,
                                                   pct_full_recovered_01 = NA_real_,
                                                   pct_full_recovered_005 = NA_real_,
                                                   full_vs_train_spearman = NA_real_,
                                                   full_vs_train_jaccard_01 = NA_real_))
    loo <- fread(loo_file)
    loo_degs_01  <- loo[padj < 0.1, gene]
    loo_degs_005 <- loo[padj < 0.05, gene]

    # Recovery
    pct_01  <- length(intersect(full_degs_01, loo_degs_01)) / max(1, length(full_degs_01)) * 100
    pct_005 <- length(intersect(full_degs_005, loo_degs_005)) / max(1, length(full_degs_005)) * 100

    # Jaccard vs full
    jac_01 <- length(intersect(full_degs_01, loo_degs_01)) /
      max(1, length(union(full_degs_01, loo_degs_01)))

    # Spearman of logFC (training vs full)
    m <- merge(full[, .(gene, full_logFC = logFC)],
               loo[, .(gene, loo_logFC = logFC)], by = "gene")
    rho <- cor(m$full_logFC, m$loo_logFC, method = "spearman",
               use = "pairwise.complete.obs")

    data.table(
      held_out_cohort         = study,
      pct_full_recovered_01   = round(pct_01, 1),
      pct_full_recovered_005  = round(pct_005, 1),
      full_vs_train_spearman  = round(rho, 4),
      full_vs_train_jaccard_01 = round(jac_01, 4)
    )
  }))

  summary_dt <- merge(summary_dt, stability_cols, by = "held_out_cohort", all.x = TRUE)
} else {
  cat("WARNING: Full dream results not found at", full_file, "\n")
  cat("Skipping stability metrics.\n")
}

# ===========================================================================
# 3. Print summary
# ===========================================================================
cat("\n===== LOO-CV v2 Summary =====\n")
print(summary_dt, digits = 4)

# Means (excluding NA)
cat("\n── Mean metrics across folds ──\n")
num_cols <- setdiff(names(summary_dt),
                    c("held_out_cohort", "n_training_cohorts"))
for (col in num_cols) {
  vals <- summary_dt[[col]]
  if (is.numeric(vals) && any(!is.na(vals))) {
    cat(sprintf("  %-30s: %.4f\n", col, mean(vals, na.rm = TRUE)))
  }
}

# ===========================================================================
# 4. Per-gene stability across folds
# ===========================================================================
cat("\n── Per-gene stability across folds ──\n")

all_loo <- rbindlist(lapply(names(metrics_list), function(study) {
  f <- file.path(LOO_DIR, paste0("dream_loo_v2_", study, ".csv"))
  if (!file.exists(f)) return(data.table())
  dt <- fread(f)
  dt[, held_out := study]
  dt
}), fill = TRUE)

if (nrow(all_loo) > 0) {
  gene_stability <- all_loo[, .(
    n_folds_tested  = .N,
    n_folds_sig_01  = sum(padj < 0.1),
    n_folds_sig_005 = sum(padj < 0.05),
    mean_logFC      = mean(logFC),
    sd_logFC        = if (.N > 1) sd(logFC) else NA_real_,
    cv_logFC        = if (.N > 1 && mean(logFC) != 0) sd(logFC) / abs(mean(logFC)) else NA_real_,
    min_padj        = min(padj),
    max_padj        = max(padj)
  ), by = gene]

  # Merge with full-model if available
  if (exists("full")) {
    gene_stability <- merge(
      gene_stability,
      full[, .(gene, full_padj = padj, full_logFC = logFC)],
      by = "gene", all = TRUE
    )
    gene_stability[, full_sig_01 := full_padj < 0.1]
    gene_stability[, full_sig_005 := full_padj < 0.05]

    # Robustness categories
    n_folds <- length(metrics_list)
    gene_stability[, robustness := fifelse(
      is.na(n_folds_tested) | n_folds_tested == 0, "Untested",
      fifelse(n_folds_sig_01 == n_folds_tested, "Robust",
              fifelse(n_folds_sig_01 >= n_folds_tested - 1, "Stable",
                      fifelse(n_folds_sig_01 >= ceiling(n_folds_tested / 2),
                              "Moderate", "Fragile")))
    )]

    cat(sprintf("Full-model DEGs (padj<0.1): %d\n",
                sum(gene_stability$full_sig_01, na.rm = TRUE)))
    cat(sprintf("  Robust (sig in all folds): %d\n",
                sum(gene_stability$full_sig_01 & gene_stability$robustness == "Robust",
                    na.rm = TRUE)))
    cat(sprintf("  Stable (sig in n-1+ folds): %d\n",
                sum(gene_stability$full_sig_01 &
                      gene_stability$robustness %in% c("Robust", "Stable"), na.rm = TRUE)))
    cat(sprintf("  Fragile (sig in <50%% folds): %d\n",
                sum(gene_stability$full_sig_01 & gene_stability$robustness == "Fragile",
                    na.rm = TRUE)))
  }

  fwrite(gene_stability, file.path(LOO_DIR, "loo_cv_per_gene_stability.csv"))
  cat("\nSaved: loo_cv_per_gene_stability.csv\n")
}

# ===========================================================================
# 5. Save summary
# ===========================================================================
fwrite(summary_dt, file.path(LOO_DIR, "loo_cv_summary.csv"))
cat("\nSaved: loo_cv_summary.csv\n")
cat("Done!\n")
