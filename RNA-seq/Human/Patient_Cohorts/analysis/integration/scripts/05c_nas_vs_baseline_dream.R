#!/usr/bin/env Rscript
# 05c_nas_vs_baseline_dream.R
# ---------------------------------------------------------------------------
# Two NAS trajectory dream mega-analyses, each using a fixed baseline (rather
# than one-vs-rest).
#
#   1. NAS_i vs NAS_0   -- baseline = samples with nas_score == 0 (no NASH
#      features per Kleiner 2005).
#   2. NAS_i vs Healthy -- baseline = samples with diagnosis_harmonized == "Control"
#      (non-MASLD), regardless of NAS score.
#
# Mirrors the structure of 60_one_vs_rest_stage_dream.R but with a fixed
# baseline instead of "rest". Only cohorts that contain BOTH the target NAS
# level AND the baseline samples participate in a given contrast. Random
# intercept per dataset.
#
# Outputs (same schema as one_vs_rest_nas_dream.csv):
#   results/staging_classifier/nas_vs_nas0_dream.csv
#   results/staging_classifier/nas_vs_healthy_dream.csv
#
# Usage: micromamba run -n rnaseq Rscript 05c_nas_vs_baseline_dream.R
# SLURM: cpu, 8 CPUs, 64GB RAM, 48h (see run_nas_baselines.sbatch)
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

cat("=== 05c: NAS vs Fixed-Baseline Dream ===\n")
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
# Helper: dream for a single NAS level vs a fixed baseline pool
# ============================================================
run_baseline_dream <- function(dge_full, meta_target, meta_baseline,
                               meta_matched_full, stage_label, param) {
  # Combine target + baseline into a single data.table with a 2-level factor
  meta_target[, arm := "target"]
  meta_baseline[, arm := "baseline"]
  meta_combined <- rbindlist(list(meta_target, meta_baseline), fill = TRUE)

  # Only cohorts that contain BOTH arms — needed for random intercept
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

  # Subset DGE
  idx <- colnames(dge_full) %in% meta_combined$sample_id
  dge_sub <- dge_full[, idx]

  # Reorder metadata to match DGE columns
  meta_ord <- meta_combined[match(colnames(dge_sub), meta_combined$sample_id)]

  # Pull inferred sex
  sex_vals <- meta_matched_full$inferred_sex[
    match(colnames(dge_sub), meta_matched_full$sample_id)]

  info <- data.frame(
    arm          = factor(meta_ord$arm, levels = c("baseline", "target")),
    dataset      = factor(meta_ord$dataset),
    inferred_sex = factor(sex_vals),
    row.names    = colnames(dge_sub),
    stringsAsFactors = FALSE
  )

  # Remove samples with missing sex
  valid <- !is.na(info$inferred_sex)
  if (sum(valid) < nrow(info)) {
    cat("  Removing", sum(!valid), "samples with missing sex\n")
    dge_sub <- dge_sub[, valid]
    info <- info[valid, , drop = FALSE]
  }
  info$dataset <- droplevels(info$dataset)

  # Re-check arm counts after sex filter
  if (sum(info$arm == "target") < 5 || sum(info$arm == "baseline") < 5) {
    cat("  SKIPPING after sex filter: insufficient samples\n")
    return(NULL)
  }
  if (nlevels(info$dataset) < 2) {
    cat("  SKIPPING: <2 datasets after sex filter\n")
    return(NULL)
  }

  # Filter + normalize
  keep <- filterByExpr(dge_sub, group = info$arm)
  dge_sub <- dge_sub[keep, , keep.lib.sizes = FALSE]
  dge_sub <- calcNormFactors(dge_sub, method = "TMM")

  # Sex may be singular within an arm/dataset; drop if only 1 level
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
# CONTRAST 1: NAS_i vs NAS_0
# ============================================================
cat("=== Contrast 1: NAS_i vs NAS_0 ===\n\n")

# NAS-scored cohorts (per CLAUDE.md + metadata audit)
NAS_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")

meta_nas <- meta_unified[
  dataset %in% NAS_DATASETS &
  !is.na(nas_score) &
  sample_id %in% pass_samples
]
meta_nas[, nas_group := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]
cat("NAS-annotated QC-passing samples:", nrow(meta_nas), "\n")
print(meta_nas[, .N, by = .(dataset, nas_group)][order(dataset, nas_group)])
cat("\n")

out1 <- file.path(OUTDIR, "nas_vs_nas0_dream.csv")
if (file.exists(out1) && file.size(out1) > 1000 && nchar(Sys.getenv("FORCE_REDO", "")) == 0) {
  cat("Already exists, skipping:", out1, "\n\n")
} else {
  baseline0 <- meta_nas[nas_group == 0]
  if (nrow(baseline0) < 5) {
    stop("Too few NAS=0 baseline samples: ", nrow(baseline0))
  }
  cat("Baseline (NAS=0) samples:", nrow(baseline0), "\n")
  cat("Baseline cohorts:", paste(unique(baseline0$dataset), collapse = ", "), "\n\n")

  all_res <- list()
  target_levels <- sort(setdiff(unique(meta_nas$nas_group), 0L))
  for (lvl in target_levels) {
    label <- paste0("NAS", lvl, "_vs_NAS0")
    cat("Running:", label, "\n")
    target <- meta_nas[nas_group == lvl]

    result <- tryCatch(
      run_baseline_dream(dge, copy(target), copy(baseline0),
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
    warning("No NAS_vs_NAS0 contrasts produced results")
  } else {
    res_dt <- rbindlist(all_res, fill = TRUE)
    setnames(res_dt, "adj.P.Val", "padj", skip_absent = TRUE)
    fwrite(res_dt, out1)
    cat("Wrote:", out1, "rows:", nrow(res_dt), "contrasts:", length(all_res), "\n\n")
  }
}

# ============================================================
# CONTRAST 2: NAS_i vs Healthy
# ============================================================
cat("=== Contrast 2: NAS_i vs Healthy ===\n\n")

# Healthy = diagnosis_harmonized == "Control". Use cohorts that have BOTH
# healthy controls AND NAS-scored disease samples. Dataset-level coverage:
#   GSE130970: 4 control samples + NAS-scored disease
#   GSE135251: 10 control + NAS-scored disease
#   GSE162694: 31 control + NAS-scored disease
# (GSE174478 and GSE193066 lack healthy controls.)
HEALTHY_DATASETS <- c("GSE130970", "GSE135251", "GSE162694")

meta_healthy <- meta_unified[
  dataset %in% HEALTHY_DATASETS &
  diagnosis_harmonized == "Control" &
  sample_id %in% pass_samples
]
meta_target_all <- meta_unified[
  dataset %in% HEALTHY_DATASETS &
  !is.na(nas_score) &
  sample_id %in% pass_samples
]
meta_target_all[, nas_group := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]

cat("Healthy (control) samples:", nrow(meta_healthy), "\n")
print(meta_healthy[, .N, by = dataset])
cat("NAS-scored samples in these cohorts:", nrow(meta_target_all), "\n\n")

out2 <- file.path(OUTDIR, "nas_vs_healthy_dream.csv")
if (file.exists(out2) && file.size(out2) > 1000 && nchar(Sys.getenv("FORCE_REDO", "")) == 0) {
  cat("Already exists, skipping:", out2, "\n\n")
} else {
  all_res <- list()
  target_levels <- sort(unique(meta_target_all$nas_group))
  for (lvl in target_levels) {
    label <- paste0("NAS", lvl, "_vs_Healthy")
    cat("Running:", label, "\n")
    target <- meta_target_all[nas_group == lvl]

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
    warning("No NAS_vs_Healthy contrasts produced results")
  } else {
    res_dt <- rbindlist(all_res, fill = TRUE)
    setnames(res_dt, "adj.P.Val", "padj", skip_absent = TRUE)
    fwrite(res_dt, out2)
    cat("Wrote:", out2, "rows:", nrow(res_dt), "contrasts:", length(all_res), "\n\n")
  }
}

cat("=== 05c_nas_vs_baseline_dream.R completed:", as.character(Sys.time()), "===\n")
