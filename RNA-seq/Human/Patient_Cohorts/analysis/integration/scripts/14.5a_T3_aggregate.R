#!/usr/bin/env Rscript
# 14.5a_T3_aggregate.R
# Aggregate per-chunk permutation results -> null_distribution.csv + empirical_pvalue.csv
# Run AFTER all array tasks finish.

suppressPackageStartupMessages({ library(data.table) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_cohort_permutation")

chunk_files <- list.files(OUT, pattern = "^chunk_\\d+_null.csv$", full.names = TRUE)
cat("Found", length(chunk_files), "chunk files\n")
if (length(chunk_files) == 0) stop("No chunk files found; nothing to aggregate.")

null_dt <- rbindlist(lapply(chunk_files, fread), use.names = TRUE, fill = TRUE)
setorder(null_dt, perm_id)
cat("Total perms aggregated:", nrow(null_dt), "\n")
fwrite(null_dt, file.path(OUT, "null_distribution.csv"))

# --- Observed value from canonical sex_interaction_dream.csv ---
obs_int <- fread(file.path(RDIR, "sex_interaction_dream.csv"))
n_obs_05 <- sum(obs_int$padj < 0.05, na.rm = TRUE)
n_obs_10 <- sum(obs_int$padj < 0.10, na.rm = TRUE)
cat("Observed n_sex_dimorphic (padj<0.05):", n_obs_05, "\n")
cat("Observed n_sex_dimorphic (padj<0.10):", n_obs_10, "\n")

null05 <- null_dt[!is.na(n_sig_pAdj05), n_sig_pAdj05]
null10 <- null_dt[!is.na(n_sig_pAdj10), n_sig_pAdj10]

emp_p_05 <- (sum(null05 >= n_obs_05) + 1L) / (length(null05) + 1L)
emp_p_10 <- (sum(null10 >= n_obs_10) + 1L) / (length(null10) + 1L)
cat("Empirical p (padj<0.05):", emp_p_05, "\n")
cat("Empirical p (padj<0.10):", emp_p_10, "\n")

emp <- data.table(
  threshold      = c("padj<0.05", "padj<0.10"),
  observed       = c(n_obs_05, n_obs_10),
  null_mean      = c(mean(null05), mean(null10)),
  null_median    = c(median(null05), median(null10)),
  null_sd        = c(sd(null05), sd(null10)),
  null_max       = c(max(null05), max(null10)),
  n_perms        = c(length(null05), length(null10)),
  empirical_pvalue = c(emp_p_05, emp_p_10)
)
fwrite(emp, file.path(OUT, "empirical_pvalue.csv"))
cat("Saved null_distribution.csv and empirical_pvalue.csv\n")
print(emp)
