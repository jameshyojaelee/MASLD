#!/usr/bin/env Rscript
# 14.4e_aggregate_bootstrap_T4.R
# ---------------------------------------------------------------------------
# Aggregate per-task bootstrap outputs from 14.4e_sex_bootstrap_T4.R into:
#   interaction_lfc_ci.csv
#     gene, n_boot_observed, boot_mean_logFC, boot_sd_logFC,
#     ci_lo_2.5, ci_hi_97.5, frac_padj_below_0.05
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ITER_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_bootstrap/iterations"
OUT_DIR  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_bootstrap"

files <- list.files(ITER_DIR, pattern = "^boot_int_t\\d+\\.csv$", full.names = TRUE)
cat("Found", length(files), "boot task files\n")
if (length(files) == 0) stop("No boot_int task files in ", ITER_DIR)

# read + bind: each task file has gene + logFC_b#### + padj_b#### columns
dts <- lapply(files, fread)
gene_ref <- dts[[1]]$gene
stopifnot(all(sapply(dts, function(d) identical(d$gene, gene_ref))))

# wide matrices
logfc_all <- do.call(cbind, lapply(dts, function(d) as.matrix(d[, grep("^logFC_b", names(d)), with = FALSE])))
padj_all  <- do.call(cbind, lapply(dts, function(d) as.matrix(d[,  grep("^padj_b",  names(d)), with = FALSE])))

cat("Total bootstrap replicates assembled:", ncol(logfc_all), "\n")

# per-gene summary
boot_mean   <- rowMeans(logfc_all, na.rm = TRUE)
boot_sd     <- apply(logfc_all, 1, sd, na.rm = TRUE)
ci_lo       <- apply(logfc_all, 1, quantile, probs = 0.025, na.rm = TRUE, names = FALSE)
ci_hi       <- apply(logfc_all, 1, quantile, probs = 0.975, na.rm = TRUE, names = FALSE)
n_obs       <- rowSums(!is.na(logfc_all))
frac_sig    <- rowMeans(padj_all < 0.05, na.rm = TRUE)

out <- data.table(
  gene                  = gene_ref,
  n_boot_observed       = n_obs,
  boot_mean_logFC       = boot_mean,
  boot_sd_logFC         = boot_sd,
  ci_lo_2_5             = ci_lo,
  ci_hi_97_5            = ci_hi,
  frac_padj_below_0_05  = frac_sig
)

fwrite(out, file.path(OUT_DIR, "interaction_lfc_ci.csv"))
cat("Wrote interaction_lfc_ci.csv (", nrow(out), " genes, B =", ncol(logfc_all), ")\n", sep = "")

# also a brief overall summary
summary_dt <- data.table(
  B_total                = ncol(logfc_all),
  n_tasks                = length(files),
  n_genes                = nrow(out),
  n_stable_sig_ge90      = sum(frac_sig >= 0.90, na.rm = TRUE),
  n_stable_sig_ge75      = sum(frac_sig >= 0.75, na.rm = TRUE),
  n_stable_sig_ge50      = sum(frac_sig >= 0.50, na.rm = TRUE),
  n_ci_excludes_zero     = sum(ci_lo > 0 | ci_hi < 0, na.rm = TRUE),
  mean_boot_sd_logFC     = mean(boot_sd, na.rm = TRUE),
  median_boot_sd_logFC   = median(boot_sd, na.rm = TRUE)
)
fwrite(summary_dt, file.path(OUT_DIR, "summary.csv"))
print(summary_dt)
cat("Done (T4 aggregation).\n")
