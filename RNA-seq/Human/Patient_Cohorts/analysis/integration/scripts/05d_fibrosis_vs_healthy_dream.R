#!/usr/bin/env Rscript
# 05d_fibrosis_vs_healthy_dream.R
# ---------------------------------------------------------------------------
# Fibrosis F_i vs Healthy dream mega-analysis — mirrors 05c (NAS_i vs Healthy).
#
# Target: samples with fibrosis_stage in {1, 2, 3, 4}
# Baseline: samples with diagnosis_harmonized == "Control" (Healthy)
#
# Cohorts participating = those that contain BOTH a given F_i AND Healthy
# controls AND at least 3 samples in each arm. Per metadata audit, this
# intersects to: GSE130970, GSE135251, GSE162694, GSE213621.
#
# Random intercept per dataset. Covariate: inferred_sex (dropped if singular).
#
# Output (schema identical to one_vs_rest_fibrosis_dream.csv and
# nas_vs_healthy_dream.csv):
#   results/staging_classifier/fibrosis_vs_healthy_dream.csv
#
# Usage: micromamba run -n rnaseq Rscript 05d_fibrosis_vs_healthy_dream.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h (see run_fibrosis_vs_healthy.sbatch)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 05d: Fibrosis vs Healthy Dream ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Parallel backend ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load data ---
cat("Loading merged counts...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))

cat("Loading matched metadata...\n")
meta_matched_file <- file.path(RDIR, "integration/meta_matched.rds")
if (!file.exists(meta_matched_file)) {
  stop("meta_matched.rds not found at: ", meta_matched_file)
}
meta_matched <- readRDS(meta_matched_file)

meta_unified <- fread(file.path(INT, "metadata/unified_metadata.csv"))
cat("Unified metadata:", nrow(meta_unified), "samples\n")

qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]
cat("QC-passing samples:", length(pass_samples), "\n\n")

# ============================================================
# Helper: dream for a single target arm vs a fixed baseline pool
# (identical structure to 05c's run_baseline_dream)
# ============================================================
run_baseline_dream <- function(dge_full, meta_target, meta_baseline,
                               meta_matched_full, stage_label, param) {
  meta_target[, arm := "target"]
  meta_baseline[, arm := "baseline"]
  meta_combined <- rbindlist(list(meta_target, meta_baseline), fill = TRUE)

  # Only cohorts that contain BOTH arms (>=3 samples each)
  per_dataset <- meta_combined[, .N, by = .(dataset, arm)]
  has_both <- dcast(per_dataset, dataset ~ arm, value.var = "N", fill = 0L)
  keep_ds <- has_both[target >= 3 & baseline >= 3, dataset]
  meta_combined <- meta_combined[dataset %in% keep_ds]

  n_target <- sum(meta_combined$arm == "target")
  n_base   <- sum(meta_combined$arm == "baseline")
  cat("  ", stage_label, ": n_target =", n_target, ", n_baseline =", n_base,
      "across", length(keep_ds), "cohort(s)\n")

  if (n_target < 5 || length(keep_ds) < 2) {
    cat("  SKIPPING: insufficient samples or <2 cohorts\n")
    return(NULL)
  }

  idx <- colnames(dge_full) %in% meta_combined$sample_id
  dge_sub <- dge_full[, idx]

  meta_ord <- meta_combined[match(colnames(dge_sub), meta_combined$sample_id)]

  sex_vals <- meta_matched_full$inferred_sex[
    match(colnames(dge_sub), meta_matched_full$sample_id)]

  info <- data.frame(
    arm          = factor(meta_ord$arm, levels = c("baseline", "target")),
    dataset      = factor(meta_ord$dataset),
    inferred_sex = factor(sex_vals),
    row.names    = colnames(dge_sub),
    stringsAsFactors = FALSE
  )

  valid <- !is.na(info$inferred_sex)
  if (sum(valid) < nrow(info)) {
    cat("  Removing", sum(!valid), "samples with missing sex\n")
    dge_sub <- dge_sub[, valid]
    info <- info[valid, , drop = FALSE]
  }
  info$dataset <- droplevels(info$dataset)

  if (sum(info$arm == "target") < 5 || sum(info$arm == "baseline") < 5) {
    cat("  SKIPPING after sex filter: insufficient samples\n")
    return(NULL)
  }
  if (nlevels(info$dataset) < 2) {
    cat("  SKIPPING: <2 datasets after sex filter\n")
    return(NULL)
  }

  keep <- filterByExpr(dge_sub, group = info$arm)
  dge_sub <- dge_sub[keep, , keep.lib.sizes = FALSE]
  dge_sub <- calcNormFactors(dge_sub, method = "TMM")

  form <- if (nlevels(info$inferred_sex) > 1) {
    ~ arm + inferred_sex + (1 | dataset)
  } else {
    ~ arm + (1 | dataset)
  }

  v <- voomWithDreamWeights(dge_sub, form, info, BPPARAM = param)
  fit <- dream(v, form, info, BPPARAM = param)

  coef_name <- "armtarget"
  if (!(coef_name %in% colnames(coef(fit)))) {
    cat("  WARNING: coefficient", coef_name, "not found\n")
    return(NULL)
  }

  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$stage_label <- stage_label

  n_sig <- sum(tt$adj.P.Val < 0.1, na.rm = TRUE)
  n_sig_strict <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
  cat("  DEGs: padj<0.1 =", n_sig, ", padj<0.05 & |LFC|>0.5 =", n_sig_strict, "\n")

  return(as.data.table(tt))
}

# ============================================================
# CONTRAST: F_i vs Healthy (i = 1..4)
# ============================================================
cat("=== Contrast: F_i vs Healthy ===\n\n")

# Per metadata audit: cohorts with BOTH healthy controls (Control) AND
# individual fibrosis stages (F0..F4).
#   GSE130970:  4 Control, F1=28/F2=9/F3=14/F4=2
#   GSE135251: 10 Control, F1=47/F2=54/F3=54/F4=14
#   GSE162694: 31 Control, F1=30/F2=27/F3=8/F4=12
#   GSE213621: 67 Control, F1=96/F2=106/F3=90  (no F4)
# GSE126848 has controls but no fibrosis staging -> excluded.
# GSE174478/GSE193066/GSE240729 have fibrosis stages but no Control arm -> excluded.
# PRJNA512027 permanently removed from pipeline 2026-05-15 (L0/S0 library batch
# perfectly confounded with disease).
FIB_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE213621")

meta_healthy <- meta_unified[
  dataset %in% FIB_DATASETS &
  diagnosis_harmonized == "Control" &
  sample_id %in% pass_samples
]
meta_target_all <- meta_unified[
  dataset %in% FIB_DATASETS &
  !is.na(fibrosis_stage) &
  as.integer(fibrosis_stage) %in% 1:4 &
  sample_id %in% pass_samples
]
meta_target_all[, fib_stage := as.integer(fibrosis_stage)]

cat("Healthy (control) samples:", nrow(meta_healthy), "\n")
print(meta_healthy[, .N, by = dataset])
cat("\nFibrosis-staged samples (F1..F4) in these cohorts:", nrow(meta_target_all), "\n")
print(meta_target_all[, .N, by = .(dataset, fib_stage)][order(dataset, fib_stage)])
cat("\n")

out <- file.path(OUTDIR, "fibrosis_vs_healthy_dream.csv")
if (file.exists(out) && file.size(out) > 1000 && nchar(Sys.getenv("FORCE_REDO", "")) == 0) {
  cat("Already exists, skipping:", out, "\n\n")
} else {
  all_res <- list()
  target_levels <- sort(unique(meta_target_all$fib_stage))
  for (lvl in target_levels) {
    label <- paste0("F", lvl, "_vs_Healthy")
    cat("Running:", label, "\n")
    target <- meta_target_all[fib_stage == lvl]

    result <- tryCatch(
      run_baseline_dream(dge, copy(target), copy(meta_healthy),
                         meta_matched, label, param),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    if (!is.null(result)) all_res[[as.character(lvl)]] <- result
    cat("\n")
  }

  if (length(all_res) == 0) {
    warning("No Fibrosis_vs_Healthy contrasts produced results")
  } else {
    res_dt <- rbindlist(all_res, fill = TRUE)
    setnames(res_dt, "adj.P.Val", "padj", skip_absent = TRUE)
    fwrite(res_dt, out)
    cat("Wrote:", out, "rows:", nrow(res_dt), "contrasts:", length(all_res), "\n\n")
  }
}

cat("=== 05d_fibrosis_vs_healthy_dream.R completed:", as.character(Sys.time()), "===\n")
