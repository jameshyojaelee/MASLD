#!/usr/bin/env Rscript
# aggregate_power_multimethod.R
# ===========================================================================
# Aggregate the 48 per-grid-cell power outputs into:
#   1. power_summary.csv     — per (grid cell x method): mean +/- sd over reps of
#        power / fdr / auroc on BOTH the own and common tested universes, plus
#        the mean tested-gene counts and mean per-rep runtime.
#   2. fdr_control_check.csv — per (grid cell x method): flags cells where the
#        mean observed FDR exceeds the 0.05 nominal threshold (anti-conservative
#        cells), for both own and common universes.
#
# Reads:  results/integration/multimethod_validation/power/power_grid_*.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR     <- file.path(PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
POWER_DIR <- file.path(RDIR, "multimethod_validation/power")

FDR_NOMINAL <- 0.05

# --- Read all grid-cell files -----------------------------------------------
# Match both the default dream+deseq2 files (power_grid_{N}.csv) and any
# single-method suffixed run (e.g. power_grid_metafor_{N}.csv) so the 3-way
# merge happens automatically. Excludes power_summary / fdr_control_check.
files <- list.files(POWER_DIR, pattern = "^power_grid.*_\\d+\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No power_grid_*.csv files in ", POWER_DIR)
cat("Reading", length(files), "grid-cell files from", POWER_DIR, "\n")

dt <- rbindlist(lapply(files, fread), fill = TRUE)
if (nrow(dt) == 0) stop("All power_grid_*.csv files were empty — nothing to aggregate")

cat("Grid cells with data:", uniqueN(dt$grid_task), "| methods:",
    paste(sort(unique(dt$method)), collapse = ", "), "\n")

GRID_KEYS <- c("grid_task", "n_per_group", "true_lfc", "tau2", "method")

# helper: mean / sd ignoring NA, with NA-safe sd (returns NA for <2 obs)
msd <- function(x) {
  x <- x[is.finite(x)]
  list(mean = if (length(x) > 0) mean(x) else NA_real_,
       sd   = if (length(x) > 1) stats::sd(x) else NA_real_)
}

# ===========================================================================
# 1. Power summary (mean +/- sd over reps)
# ===========================================================================
summ <- dt[, {
  po  <- msd(power_own);    pc  <- msd(power_common)
  fo  <- msd(fdr_own);      fc  <- msd(fdr_common)
  ao  <- msd(auroc_own);    ac  <- msd(auroc_common)
  .(
    n_reps          = .N,
    n_truth_de_mean = mean(n_truth_de, na.rm = TRUE),
    n_tested_mean        = mean(n_tested, na.rm = TRUE),
    n_tested_common_mean = mean(n_tested_common, na.rm = TRUE),
    power_own_mean  = po$mean, power_own_sd  = po$sd,
    fdr_own_mean    = fo$mean, fdr_own_sd    = fo$sd,
    auroc_own_mean  = ao$mean, auroc_own_sd  = ao$sd,
    power_common_mean = pc$mean, power_common_sd = pc$sd,
    fdr_common_mean   = fc$mean, fdr_common_sd   = fc$sd,
    auroc_common_mean = ac$mean, auroc_common_sd = ac$sd,
    runtime_min_mean  = mean(runtime_min, na.rm = TRUE)
  )
}, by = GRID_KEYS]

# round numeric summaries for readability (keep counts as-is)
num_round <- setdiff(names(summ),
  c(GRID_KEYS, "n_reps", "n_truth_de_mean",
    "n_tested_mean", "n_tested_common_mean", "runtime_min_mean"))
summ[, (num_round) := lapply(.SD, function(x) round(x, 4)), .SDcols = num_round]
summ[, `:=`(n_truth_de_mean = round(n_truth_de_mean, 1),
            n_tested_mean = round(n_tested_mean, 1),
            n_tested_common_mean = round(n_tested_common_mean, 1),
            runtime_min_mean = round(runtime_min_mean, 3))]
setorder(summ, grid_task, method)
fwrite(summ, file.path(POWER_DIR, "power_summary.csv"))
cat("Wrote power_summary.csv:", nrow(summ), "grid-cell x method rows\n")

# ===========================================================================
# 2. FDR control check — flag anti-conservative cells (mean observed FDR > 0.05)
# ===========================================================================
fdr_check <- summ[, .(
  grid_task, n_per_group, true_lfc, tau2, method,
  fdr_own_mean, fdr_common_mean,
  fdr_own_exceeds    = is.finite(fdr_own_mean)    & fdr_own_mean    > FDR_NOMINAL,
  fdr_common_exceeds = is.finite(fdr_common_mean) & fdr_common_mean > FDR_NOMINAL
)]
setorder(fdr_check, -fdr_own_mean)
fwrite(fdr_check, file.path(POWER_DIR, "fdr_control_check.csv"))

n_own_flag    <- sum(fdr_check$fdr_own_exceeds,    na.rm = TRUE)
n_common_flag <- sum(fdr_check$fdr_common_exceeds, na.rm = TRUE)
cat("Wrote fdr_control_check.csv:", nrow(fdr_check), "rows | anti-conservative cells: ",
    n_own_flag, "(own) /", n_common_flag, "(common)\n")
cat("Aggregation complete.\n")
