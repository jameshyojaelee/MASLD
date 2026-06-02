#!/usr/bin/env Rscript
# lfc_grid_sweep.R
# ---------------------------------------------------------------------------
# B3 — LFC × padj grid sweep with 5-fold LOO recovery (empirical plateau).
#
# For every (lfc, padj) combination in:
#     lfc  ∈ {0.10, 0.15, 0.20, ..., 1.00}  (19 values)
#     padj ∈ {0.01, 0.05, 0.10}             (3  values)
# we recompute:
#   * n_DEG = #{genes : padj_full < padj & |logFC_full| > lfc}
#   * per-fold recovery (5 folds, holding out each of the 5 mega cohorts):
#         recovery_k = |full_DEG ∩ loo_DEG_k| / |full_DEG|
#   * mean_LOO_recovery = mean(recovery_k)
#
# Reuses pre-computed dream LOO results in
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv/
# (one dream_loo_<COHORT>.csv per mega cohort, all 5 already present).
# This avoids re-running dream() 5×60 = 300 times.
#
# Output: RNA-seq/results/audit_sensitivity/lfc_sweep/lfc_sweep_results.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

# Allow override for worktree execution
if (!dir.exists(BASE)) {
  stop("MASLD_PROJECT_ROOT does not exist: ", BASE)
}

INT_RES <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR <- file.path(INT_RES, "loo_cv")
FULL_F  <- file.path(INT_RES, "dream_results.csv")

OUT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/lfc_sweep")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("B3: LFC × padj grid sweep with 5-fold LOO recovery\n")
cat("============================================================\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("Base:", BASE, "\n")
cat("Full dream:", FULL_F, "\n")
cat("LOO dir   :", LOO_DIR, "\n\n")

stopifnot(file.exists(FULL_F))
stopifnot(dir.exists(LOO_DIR))

full <- fread(FULL_F)
cat("Full dream rows:", nrow(full), "\n")
stopifnot(all(c("gene", "logFC", "padj") %in% names(full)))

# 5 mega cohorts (those with healthy controls; matches yaml include_in_mega).
# LOO files are named dream_loo_<COHORT>.csv with these accession IDs.
MEGA_COHORTS <- c("GSE126848", "GSE130970", "GSE135251",
                  "GSE162694", "GSE213621")

loo_list <- list()
for (co in MEGA_COHORTS) {
  f <- file.path(LOO_DIR, paste0("dream_loo_", co, ".csv"))
  if (!file.exists(f)) stop("Missing LOO dream result: ", f)
  d <- fread(f)
  stopifnot(all(c("gene", "logFC", "padj") %in% names(d)))
  loo_list[[co]] <- d[, .(gene, logFC, padj)]
  cat("Loaded LOO", co, ":", nrow(d), "rows\n")
}
cat("\n")

# Grid
lfc_vals  <- seq(0.10, 1.00, by = 0.05)   # 19 values
padj_vals <- c(0.01, 0.05, 0.10)          # 3 values
cat("LFC values  (", length(lfc_vals), "):", paste(lfc_vals, collapse = ", "), "\n")
cat("padj values (", length(padj_vals), "):", paste(padj_vals, collapse = ", "), "\n")
cat("Total combinations:", length(lfc_vals) * length(padj_vals), "\n\n")

results <- list()
i <- 0
for (lfc in lfc_vals) {
  for (pa in padj_vals) {
    i <- i + 1
    full_sig <- !is.na(full$padj) & !is.na(full$logFC) &
                full$padj < pa & abs(full$logFC) > lfc
    full_degs <- full$gene[full_sig]
    n_full <- length(full_degs)

    per_fold <- numeric(length(MEGA_COHORTS))
    names(per_fold) <- MEGA_COHORTS

    for (co in MEGA_COHORTS) {
      d <- loo_list[[co]]
      loo_sig <- !is.na(d$padj) & !is.na(d$logFC) &
                 d$padj < pa & abs(d$logFC) > lfc
      loo_degs <- d$gene[loo_sig]
      per_fold[co] <- if (n_full == 0) NA_real_ else
                       length(intersect(full_degs, loo_degs)) / n_full
    }

    mean_rec <- mean(per_fold, na.rm = TRUE)
    sd_rec   <- sd(per_fold,   na.rm = TRUE)

    results[[i]] <- data.table(
      lfc                  = lfc,
      padj                 = pa,
      n_DEG                = n_full,
      mean_LOO_recovery    = mean_rec,
      sd_LOO_recovery      = sd_rec,
      recovery_GSE126848   = per_fold["GSE126848"],
      recovery_GSE130970   = per_fold["GSE130970"],
      recovery_GSE135251   = per_fold["GSE135251"],
      recovery_GSE162694   = per_fold["GSE162694"],
      recovery_GSE213621   = per_fold["GSE213621"]
    )
  }
}

res <- rbindlist(results)
cat("Grid computed. Head:\n")
print(head(res, 12))
cat("\nTail:\n")
print(tail(res, 6))

out_f <- file.path(OUT_DIR, "lfc_sweep_results.csv")
fwrite(res, out_f)
cat("\nSaved:", out_f, "\n")

# Quick plateau detection: derivative of mean recovery vs LFC at padj=0.05
sub <- res[padj == 0.05][order(lfc)]
sub[, drec := c(NA, diff(mean_LOO_recovery))]
cat("\nPlateau diagnostic at padj=0.05:\n")
print(sub[, .(lfc, n_DEG, mean_LOO_recovery = round(mean_LOO_recovery, 4),
              drec = round(drec, 4))])

cat("\nDone at", as.character(Sys.time()), "\n")
