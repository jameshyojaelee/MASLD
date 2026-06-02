#!/usr/bin/env Rscript
# aggregate_pillar_C.R
# ---------------------------------------------------------------------------
# Aggregate Pillar C (within-cohort permutation null) per-iter outputs.
# Computes:
#   * Distribution of permuted DEG counts (mean / median / 95th percentile)
#   * Empirical FDR = median permuted count / canonical DEG count (4,370, kallisto canonical)
#   * Per-gene permutation z-score = (observed t - mean(perm t)) / sd(perm t)
#   * Per-gene empirical p-value = P(|perm t| >= |observed t|)
#
# Out: audit_sensitivity/pillar_C_empirical_null.csv      (per-gene)
#      audit_sensitivity/pillar_C_perm_count_dist.csv     (per-iter)
#      audit_sensitivity/pillar_C_summary.csv             (1 row, headline numbers)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
PERM_DIR <- file.path(OUT_DIR, "permutation_null")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ_THR <- 0.05
LFC_THR  <- as.numeric(Sys.getenv("LFC_THR", "0.5"))
THR_TAG  <- if (LFC_THR == 0) "nolfc" else paste0("lfc", sub("^0\\.", "", as.character(LFC_THR)))
cat(sprintf("Aggregating Pillar C at padj<%.2f & |logFC|>%.2f  (tag = %s)\n",
            PADJ_THR, LFC_THR, THR_TAG))

# --- Observed (primary dream) ---
full <- fread(file.path(RDIR, "dream_results.csv"))
deg_full <- full[padj < PADJ_THR & abs(logFC) > LFC_THR, gene]
N_DEG_OBS <- length(deg_full)
cat(sprintf("Observed DEGs (canonical): %d\n", N_DEG_OBS))

perm_files <- list.files(PERM_DIR, pattern = "^perm_iter_\\d{4}\\.csv$", full.names = TRUE)
cat(sprintf("Permutation files: %d (target B=1000)\n", length(perm_files)))
if (length(perm_files) == 0) stop("No permutation files found; run dream_disease_permutation_iter.R first.")

# ---- Per-iter DEG count + per-gene t-statistic stack ----
perm_long <- rbindlist(lapply(perm_files, function(f) {
  dt <- fread(f, select = c("gene", "logFC", "t", "padj"))
  dt[, iter := as.integer(sub("^perm_iter_(\\d{4})\\.csv$", "\\1", basename(f)))]
  dt
}))

# DEG counts per permutation
count_dt <- perm_long[, .(n_deg_perm = sum(padj < PADJ_THR & abs(logFC) > LFC_THR)),
                      by = iter]
fwrite(count_dt, file.path(OUT_DIR, paste0("pillar_C_perm_count_dist_", THR_TAG, ".csv")))

n_perm <- nrow(count_dt)
median_perm <- median(count_dt$n_deg_perm)
mean_perm   <- mean(count_dt$n_deg_perm)
q95_perm    <- quantile(count_dt$n_deg_perm, 0.95)
emp_fdr     <- median_perm / N_DEG_OBS

cat(sprintf("\nPermuted DEG count: median=%.1f mean=%.1f q95=%.1f (B=%d)\n",
            median_perm, mean_perm, q95_perm, n_perm))
cat(sprintf("Empirical FDR (median permuted / observed) = %.4f\n", emp_fdr))
cat(sprintf("Pass criterion (emp FDR < 0.05): %s\n",
            ifelse(emp_fdr < 0.05, "PASS", "FAIL")))

# ---- Per-gene null distribution → z-score + empirical p ----
# Observed t-statistic comes from the full dream fit. We use absolute value to be
# two-sided. Need the observed t for each gene.
obs_t <- full[, .(gene, t_obs = t)]

# Null statistics per gene (over all permutations)
null_stats <- perm_long[, .(t_perm_mean = mean(t),
                            t_perm_sd   = sd(t),
                            n_perm_per_gene = .N), by = gene]

# Empirical two-sided p-value: P(|perm t| >= |obs t|)
perm_long_obs <- merge(perm_long, obs_t, by = "gene")
perm_long_obs[, exceed := abs(t) >= abs(t_obs)]
emp_p <- perm_long_obs[, .(perm_emp_p = mean(exceed)), by = gene]

per_gene <- merge(obs_t, null_stats, by = "gene", all.x = TRUE)
per_gene <- merge(per_gene, emp_p,    by = "gene", all.x = TRUE)
per_gene[, perm_z := round((t_obs - t_perm_mean) / t_perm_sd, 3)]
per_gene[, perm_emp_p := round(perm_emp_p, 5)]
per_gene[, t_obs := round(t_obs, 4)]
per_gene[, t_perm_mean := round(t_perm_mean, 4)]
per_gene[, t_perm_sd   := round(t_perm_sd, 4)]
fwrite(per_gene, file.path(OUT_DIR, paste0("pillar_C_empirical_null_", THR_TAG, ".csv")))
if (THR_TAG == "lfc05") fwrite(per_gene, file.path(OUT_DIR, "pillar_C_empirical_null.csv"))
cat("Saved per-gene null:", nrow(per_gene), "rows\n")

# ---- Headline summary ----
n_canonical_passing <- per_gene[gene %in% deg_full & abs(perm_z) > 3, .N]
sum_dt <- data.table(
  n_perm        = n_perm,
  n_deg_observed = N_DEG_OBS,
  perm_count_median = round(median_perm, 1),
  perm_count_mean   = round(mean_perm, 1),
  perm_count_q95    = round(as.numeric(q95_perm), 1),
  empirical_fdr     = round(emp_fdr, 4),
  pass_fdr_lt_5pct  = emp_fdr < 0.05,
  n_canonical_perm_z_gt_3 = n_canonical_passing,
  pct_canonical_perm_z_gt_3 = round(100 * n_canonical_passing / N_DEG_OBS, 2))
fwrite(sum_dt, file.path(OUT_DIR, paste0("pillar_C_summary_", THR_TAG, ".csv")))
if (THR_TAG == "lfc05") fwrite(sum_dt, file.path(OUT_DIR, "pillar_C_summary.csv"))
cat("\n=== Pillar C summary ===\n"); print(sum_dt)
cat("\nDone Pillar C aggregation.\n")
