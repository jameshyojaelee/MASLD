#!/usr/bin/env Rscript
# ===========================================================================
# aggregate_lvqw_loo_cv_C2.R
# Aggregate the LVQW LOO-CV (C2) folds. Recover each fold's training DEGs
# against the FULL C2 canonical (canonical_deg_results.csv) — mirrors
# aggregate_loo_cv_v2.R's v1-style "stability" metrics but pointed at the C2
# canonical instead of dream_results.csv.
#
# Reports (means across 5 folds):
#   mean recovery (% of full C2 DEGs recovered; padj<.05 and padj<.1)
#   full-vs-train Spearman rho, direction concordance, Jaccard
#   AUC replication (held-out), and the ROBUST gene count (full DEG sig in all folds).
#
# Output: results/integration/loo_cv_C2/loo_cv_C2_summary.csv
#         results/integration/loo_cv_C2/loo_cv_C2_per_gene_stability.csv
# ===========================================================================
suppressPackageStartupMessages({ library(data.table); library(yaml) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO  <- file.path(RDIR, "loo_cv_C2")

ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
COHORTS <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
cat("LOO cohort set (k =", length(COHORTS), "):", paste(COHORTS, collapse = ", "), "\n")

# --- per-fold metrics (held-out replication) ---
metrics_list <- list(); missing <- character(0)
for (study in COHORTS) {
  f <- file.path(LOO, paste0("loo_C2_metrics_", study, ".csv"))
  if (!file.exists(f)) { missing <- c(missing, study); next }
  metrics_list[[study]] <- fread(f)
}
if (length(missing)) cat("*** missing folds:", paste(missing, collapse = ", "), "***\n")
if (!length(metrics_list)) stop("No LVQW LOO-CV (C2) fold metrics found!")
summary_dt <- rbindlist(metrics_list, fill = TRUE)

# --- recovery vs FULL C2 canonical ---
full <- fread(file.path(RDIR, "canonical_deg_results.csv"))
full_degs_01  <- full[padj < 0.1, gene]
full_degs_005 <- full[padj < 0.05, gene]
full_degs_lfc <- full[padj < 0.05 & abs(logFC) > 0.5, gene]
full_degs_ash <- full[!is.na(lfsr) & lfsr < 0.05 & abs(shrunk_logFC) > 0.5, gene]
cat(sprintf("Full C2 canonical: %d DEGs (padj<.1), %d (padj<.05), %d raw-Tier1, %d ashr-Tier1\n",
            length(full_degs_01), length(full_degs_005), length(full_degs_lfc), length(full_degs_ash)))

stability_cols <- rbindlist(lapply(names(metrics_list), function(study) {
  loo_file <- file.path(LOO, paste0("lvqw_loo_C2_", study, ".csv"))
  if (!file.exists(loo_file)) return(data.table(held_out_cohort = study))
  loo <- fread(loo_file)
  loo_degs_01  <- loo[padj < 0.1, gene]
  loo_degs_005 <- loo[padj < 0.05, gene]
  loo_degs_ash <- loo[!is.na(lfsr) & lfsr < 0.05 & abs(shrunk_logFC) > 0.5, gene]

  pct_01  <- length(intersect(full_degs_01, loo_degs_01))   / max(1, length(full_degs_01))  * 100
  pct_005 <- length(intersect(full_degs_005, loo_degs_005)) / max(1, length(full_degs_005)) * 100
  pct_ash <- length(intersect(full_degs_ash, loo_degs_ash)) / max(1, length(full_degs_ash)) * 100
  jac_01  <- length(intersect(full_degs_01, loo_degs_01)) / max(1, length(union(full_degs_01, loo_degs_01)))

  m   <- merge(full[, .(gene, full_logFC = logFC)], loo[, .(gene, loo_logFC = logFC)], by = "gene")
  rho <- cor(m$full_logFC, m$loo_logFC, method = "spearman", use = "pairwise.complete.obs")
  # direction concordance among full DEGs (raw Tier1) recovered in this fold's logFC
  md  <- m[gene %in% full_degs_lfc]
  dir_full <- if (nrow(md) > 0) mean(sign(md$full_logFC) == sign(md$loo_logFC)) * 100 else NA_real_

  data.table(held_out_cohort = study,
             pct_full_recovered_01 = round(pct_01, 1),
             pct_full_recovered_005 = round(pct_005, 1),
             pct_full_recovered_ashTier1 = round(pct_ash, 1),
             full_vs_train_spearman = round(rho, 4),
             full_vs_train_jaccard_01 = round(jac_01, 4),
             full_vs_train_direction = round(dir_full, 1))
}), fill = TRUE)

summary_dt <- merge(summary_dt, stability_cols, by = "held_out_cohort", all.x = TRUE)

cat("\n===== LVQW LOO-CV (C2) per-fold summary =====\n")
print(summary_dt[, .(held_out_cohort, n_samples_training, gene_universe_size,
                     n_degs_training_005, pct_full_recovered_005, pct_full_recovered_01,
                     full_vs_train_spearman, full_vs_train_direction, full_vs_train_jaccard_01,
                     jaccard_01, lfc_spearman, direction_concordance, auc_replication)], digits = 4)

cat("\n-- Mean metrics across folds --\n")
mean_metrics <- c("pct_full_recovered_005", "pct_full_recovered_01", "pct_full_recovered_ashTier1",
                  "full_vs_train_spearman", "full_vs_train_direction", "full_vs_train_jaccard_01",
                  "jaccard_01", "lfc_spearman", "lfc_spearman_sig", "direction_concordance",
                  "auc_replication", "fisher_or")
for (col in mean_metrics) {
  if (col %in% names(summary_dt) && any(!is.na(summary_dt[[col]])))
    cat(sprintf("  %-30s: %.4f\n", col, mean(summary_dt[[col]], na.rm = TRUE)))
}

# --- per-gene stability across folds + ROBUST count ---
all_loo <- rbindlist(lapply(names(metrics_list), function(study) {
  f <- file.path(LOO, paste0("lvqw_loo_C2_", study, ".csv"))
  if (!file.exists(f)) return(data.table())
  dt <- fread(f); dt[, held_out := study]; dt
}), fill = TRUE)

if (nrow(all_loo) > 0) {
  gene_stability <- all_loo[, .(
    n_folds_tested  = .N,
    n_folds_sig_01  = sum(padj < 0.1),
    n_folds_sig_005 = sum(padj < 0.05)
  ), by = gene]
  gene_stability <- merge(gene_stability,
                          full[, .(gene, full_padj = padj, full_logFC = logFC)],
                          by = "gene", all = TRUE)
  gene_stability[, full_sig_01 := full_padj < 0.1]
  gene_stability[, full_sig_005 := full_padj < 0.05]
  n_folds <- length(metrics_list)
  gene_stability[, robustness := fifelse(
    is.na(n_folds_tested) | n_folds_tested == 0, "Untested",
    fifelse(n_folds_sig_01 == n_folds_tested, "Robust",
            fifelse(n_folds_sig_01 >= n_folds_tested - 1, "Stable",
                    fifelse(n_folds_sig_01 >= ceiling(n_folds_tested / 2), "Moderate", "Fragile"))))]

  n_robust <- sum(gene_stability$full_sig_01 & gene_stability$robustness == "Robust", na.rm = TRUE)
  n_stable <- sum(gene_stability$full_sig_01 & gene_stability$robustness %in% c("Robust", "Stable"), na.rm = TRUE)
  cat(sprintf("\nFull-model DEGs (padj<.1): %d\n", sum(gene_stability$full_sig_01, na.rm = TRUE)))
  cat(sprintf("  Robust (sig in all %d folds): %d\n", n_folds, n_robust))
  cat(sprintf("  Stable (sig in n-1+ folds):   %d\n", n_stable))
  fwrite(gene_stability, file.path(LOO, "loo_cv_C2_per_gene_stability.csv"))
}

fwrite(summary_dt, file.path(LOO, "loo_cv_C2_summary.csv"))
cat("\nSaved: loo_cv_C2_summary.csv + loo_cv_C2_per_gene_stability.csv\n")
cat("Done!\n")
