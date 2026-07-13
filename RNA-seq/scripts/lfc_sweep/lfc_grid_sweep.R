#!/usr/bin/env Rscript
# lfc_grid_sweep.R
# ---------------------------------------------------------------------------
# B3 — LFC × significance grid sweep with 5-fold LOO recovery (empirical plateau).
#
# Two scales are produced:
#   raw    : gate  padj < {0.01,0.05,0.10}  &  |logFC|        > lfc
#   shrunk : gate  lfsr < {0.01,0.05,0.10}  &  |shrunk_logFC| > lfc
#
# The full dream model + the 5 LOO dream refits are RAW on disk. The shrunk
# scale is obtained by applying ashr on the fly to each table using
# SE (full carries an SE column; folds derive SE = |logFC / t|). Keeping both
# the full model and the folds on the SAME (dream) method makes the LOO
# recovery numerator/denominator method-consistent; the canonical
# limma-voom-qw C2 table is NOT substituted here because the on-disk LOO folds
# are dream refits (re-running limma-voom LOO is out of scope).
#
# For every (lfc, sig) combination:
#     lfc  ∈ {0.10, 0.15, ..., 1.00}        (19 values)
#     sig  ∈ {0.01, 0.05, 0.10}             (3  values)
# we recompute:
#   * n_DEG = #{genes : sig_full < sig & |eff_full| > lfc}
#   * per-fold recovery (5 folds, holding out each of the 5 mega cohorts):
#         recovery_k = |full_DEG ∩ loo_DEG_k| / |full_DEG|
#   * mean_LOO_recovery = mean(recovery_k)
#
# Reuses pre-computed dream LOO results in
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv/
# (one dream_loo_<COHORT>.csv per mega cohort, all 5 already present).
#
# Output:
#   RNA-seq/results/audit_sensitivity/lfc_sweep/lfc_sweep_results_raw.csv
#   RNA-seq/results/audit_sensitivity/lfc_sweep/lfc_sweep_results_shrunk.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ashr)
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
cat("B3: LFC × significance grid sweep with 5-fold LOO recovery\n")
cat("============================================================\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("Base:", BASE, "\n")
cat("Full dream:", FULL_F, "\n")
cat("LOO dir   :", LOO_DIR, "\n\n")

stopifnot(file.exists(FULL_F))
stopifnot(dir.exists(LOO_DIR))

# --- ashr helper: shrink (logFC, SE) -> shrunk_logFC + lfsr -----------------
shrink_dt <- function(d) {
  d <- copy(d)
  if (!"SE" %in% names(d)) {
    stopifnot("t" %in% names(d))
    d[, SE := abs(logFC / t)]
  }
  ok <- is.finite(d$logFC) & is.finite(d$SE) & d$SE > 0
  fit <- ashr::ash(d$logFC[ok], d$SE[ok], mixcompdist = "normal")
  d[, shrunk_logFC := NA_real_][, lfsr := NA_real_]
  d[ok, shrunk_logFC := ashr::get_pm(fit)]
  d[ok, lfsr := ashr::get_lfsr(fit)]
  d[]
}

full <- fread(FULL_F)
cat("Full dream rows:", nrow(full), "\n")
stopifnot(all(c("gene", "logFC", "padj") %in% names(full)))
full <- shrink_dt(full)

# 5 mega cohorts (those with healthy controls; matches yaml include_in_mega).
MEGA_COHORTS <- c("GSE126848", "GSE130970", "GSE135251",
                  "GSE162694", "GSE213621")

loo_list <- list()
for (co in MEGA_COHORTS) {
  f <- file.path(LOO_DIR, paste0("dream_loo_", co, ".csv"))
  if (!file.exists(f)) stop("Missing LOO dream result: ", f)
  d <- fread(f)
  stopifnot(all(c("gene", "logFC", "padj", "t") %in% names(d)))
  loo_list[[co]] <- shrink_dt(d[, .(gene, logFC, padj, t)])   # adds shrunk_logFC + lfsr
  cat("Loaded LOO", co, ":", nrow(d), "rows (shrunk via SE=|logFC/t|)\n")
}
cat("\n")

# Grid
lfc_vals <- seq(0.10, 1.00, by = 0.05)   # 19 values
sig_vals <- c(0.01, 0.05, 0.10)          # 3 values

SCALES <- list(
  raw    = list(eff = "logFC",        sig = "padj", sig_type = "padj", suffix = "raw"),
  shrunk = list(eff = "shrunk_logFC", sig = "lfsr", sig_type = "lfsr", suffix = "shrunk")
)

for (scale in names(SCALES)) {
  sc  <- SCALES[[scale]]
  eff <- sc$eff; sig <- sc$sig
  cat(sprintf("--- scale = %s (gate %s < sig & |%s| > lfc) ---\n", scale, sig, eff))

  results <- list(); i <- 0
  for (lfc in lfc_vals) {
    for (sg in sig_vals) {
      i <- i + 1
      full_sig <- !is.na(full[[sig]]) & !is.na(full[[eff]]) &
                  full[[sig]] < sg & abs(full[[eff]]) > lfc
      full_degs <- full$gene[full_sig]
      n_full <- length(full_degs)

      per_fold <- numeric(length(MEGA_COHORTS)); names(per_fold) <- MEGA_COHORTS
      for (co in MEGA_COHORTS) {
        d <- loo_list[[co]]
        loo_sig <- !is.na(d[[sig]]) & !is.na(d[[eff]]) &
                   d[[sig]] < sg & abs(d[[eff]]) > lfc
        loo_degs <- d$gene[loo_sig]
        per_fold[co] <- if (n_full == 0) NA_real_ else
                         length(intersect(full_degs, loo_degs)) / n_full
      }

      results[[i]] <- data.table(
        scale                = scale,
        sig_type             = sc$sig_type,
        lfc                  = lfc,
        sig_thr              = sg,
        n_DEG                = n_full,
        mean_LOO_recovery    = mean(per_fold, na.rm = TRUE),
        sd_LOO_recovery      = sd(per_fold,   na.rm = TRUE),
        recovery_GSE126848   = per_fold["GSE126848"],
        recovery_GSE130970   = per_fold["GSE130970"],
        recovery_GSE135251   = per_fold["GSE135251"],
        recovery_GSE162694   = per_fold["GSE162694"],
        recovery_GSE213621   = per_fold["GSE213621"]
      )
    }
  }
  res <- rbindlist(results)

  out_f <- file.path(OUT_DIR, sprintf("lfc_sweep_results_%s.csv", sc$suffix))
  fwrite(res, out_f)
  cat("Saved:", out_f, "\n")

  # Plateau diagnostic at sig = 0.05
  sub <- res[sig_thr == 0.05][order(lfc)]
  sub[, drec := c(NA, diff(mean_LOO_recovery))]
  cat(sprintf("Plateau diagnostic [%s] at %s<0.05:\n", scale, sc$sig_type))
  print(sub[, .(lfc, n_DEG, mean_LOO_recovery = round(mean_LOO_recovery, 4),
                drec = round(drec, 4))])
  cat("\n")
}

cat("Done at", as.character(Sys.time()), "\n")
