#!/usr/bin/env Rscript
# 14b_stage_vs_healthy_dream.R
# Per-fibrosis-stage and per-NAS-score dream mega-analyses using STRICTLY
# HEALTHY controls (condition == "Control") as reference, not all F0/NAS0 samples.
#
# Motivation: The existing 14b uses all F0/NAS0 samples as reference, which
# includes NAFL/NASH_F0 patients. This analysis uses only the 64 explicitly
# labeled healthy control samples (GSE130970 n=23, GSE135251 n=10, GSE162694 n=31)
# as reference, giving cleaner disease-vs-healthy contrasts.
#
# Produces:
#   fibrosis_stage_vs_ctrl_dream.csv   (F1_vs_Ctrl, F2_vs_Ctrl, F3_vs_Ctrl, F4_vs_Ctrl)
#   nas_stage_vs_ctrl_dream.csv        (NAS1_vs_Ctrl ... NAS7_vs_Ctrl)
#   fibrosis_stage_vs_ctrl_sample_sizes.csv
#   nas_stage_vs_ctrl_sample_sizes.csv
#
# Usage: Rscript 14b_stage_vs_healthy_dream.R
# SLURM: bigmem, 16 CPU, 200GB RAM, ~2-4h

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
SIGS <- file.path(RDIR, "disease_signatures")
dir.create(SIGS, showWarnings = FALSE, recursive = TRUE)

cat("=== 14b_stage_vs_healthy_dream: Stage DEGs vs Strictly Healthy Controls ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Setup parallel backend ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load data ---
cat("Loading merged DGE and metadata...\n")
dge      <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))
meta_mat <- readRDS(file.path(RDIR, "integration/meta_matched.rds"))
meta     <- fread(file.path(INT, "metadata/unified_metadata.csv"))
qc       <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]
cat("QC-passing samples:", length(pass_samples), "\n\n")

# Helper: run dream stage contrasts
run_stage_dream <- function(dge_full, info, form, coef_prefix, stage_levels,
                             contrast_fmt, param) {
  # Expects info to have a `grp` factor column with "Ctrl" as reference
  keep <- filterByExpr(dge_full, group = info$grp)
  dge_sub <- dge_full[keep, , keep.lib.sizes = FALSE]
  dge_sub  <- calcNormFactors(dge_sub, method = "TMM")
  cat("  Genes after expression filter:", nrow(dge_sub), "\n")

  cat("  Running voomWithDreamWeights...\n")
  v <- voomWithDreamWeights(dge_sub, form, info, BPPARAM = param)

  cat("  Running dream...\n")
  fit <- dream(v, form, info, BPPARAM = param)
  # DO NOT call eBayes() — dream() already applies moderated t-stats

  results <- list()
  for (lvl in stage_levels) {
    coef_name <- paste0(coef_prefix, lvl)
    if (!coef_name %in% colnames(coef(fit))) {
      cat("  WARNING: coefficient", coef_name, "not found — skipping\n")
      next
    }
    tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
    tt$gene     <- rownames(tt)
    tt$contrast <- sprintf(contrast_fmt, lvl)
    results[[as.character(lvl)]] <- as.data.table(tt)
    n_sig <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
    cat(sprintf("  %s: %d DEGs (padj<0.05, |LFC|>0.5)\n",
                sprintf(contrast_fmt, lvl), n_sig))
  }
  rbindlist(results, fill = TRUE)
}

# ============================================================
# ANALYSIS 1: Per-Fibrosis-Stage Dream vs Healthy Controls
# ============================================================
cat("=== ANALYSIS 1: Fibrosis Stages vs Healthy Controls ===\n")

FIB_DATASETS <- c("GSE130970", "GSE135251", "GSE162694",
                   "GSE174478", "GSE193066", "GSE240729")

meta_fib_all <- meta[
  dataset %in% FIB_DATASETS &
  sample_id %in% pass_samples
]

# Control group: condition == "Control"
# Disease groups: samples with fibrosis_stage 1-4
# Exclude: F0 non-control samples (NAFL/NASH_F0/Fibrosis_F0 — they are
#          neither strictly healthy nor a discrete fibrosis stage)
meta_ctrl_fib <- meta_fib_all[condition == "Control"]
meta_ctrl_fib[, grp := "Ctrl"]

meta_fib_dis <- meta_fib_all[
  !is.na(fibrosis_stage) &
  as.integer(fibrosis_stage) %in% 1:4 &
  condition != "Control"
]
meta_fib_dis[, grp := paste0("F", as.integer(fibrosis_stage))]

meta_fib <- rbind(meta_ctrl_fib, meta_fib_dis, fill = TRUE)
cat("Reference (Control):", nrow(meta_ctrl_fib), "\n")
print(meta_ctrl_fib[, .N, by = dataset][order(dataset)])

fib_dist <- meta_fib[, .N, by = grp][order(grp)]
cat("\nSample distribution (group):\n")
print(fib_dist)
fwrite(fib_dist, file.path(SIGS, "fibrosis_stage_vs_ctrl_sample_sizes.csv"))

# Subset DGE
idx_fib <- colnames(dge) %in% meta_fib$sample_id
dge_fib  <- dge[, idx_fib]
cat("DGE subset:", ncol(dge_fib), "samples\n")

# Build info data.frame aligned to DGE columns
idx_order    <- match(colnames(dge_fib), meta_fib$sample_id)
matched_fib  <- meta_fib[idx_order]
inferred_sex <- meta_mat$inferred_sex[match(colnames(dge_fib), meta_mat$sample_id)]

info_fib <- data.frame(
  grp          = factor(matched_fib$grp, levels = c("Ctrl", "F1", "F2", "F3", "F4")),
  dataset      = factor(matched_fib$dataset),
  inferred_sex = factor(inferred_sex),
  row.names    = colnames(dge_fib),
  stringsAsFactors = FALSE
)

cat("  grp level counts:\n")
print(table(info_fib$grp, useNA = "always"))

# Require at least 3 datasets for the random effect to be estimable
n_datasets <- length(unique(info_fib$dataset[info_fib$grp != "Ctrl"]))
cat("  Datasets contributing disease samples:", n_datasets, "\n")

form_fib <- ~ grp + inferred_sex + (1 | dataset)
cat("\nDream formula:", deparse(form_fib), "\n")

fib_results <- run_stage_dream(
  dge_fib, info_fib,
  form    = form_fib,
  coef_prefix   = "grp",
  stage_levels  = paste0("F", 1:4),
  contrast_fmt  = "%s_vs_Ctrl",
  param   = param
)
setnames(fib_results, "adj.P.Val", "padj", skip_absent = TRUE)
fwrite(fib_results, file.path(SIGS, "fibrosis_stage_vs_ctrl_dream.csv"))
cat("\nFibrosis vs control dream results saved:", nrow(fib_results), "rows\n")

# ============================================================
# ANALYSIS 2: Per-NAS-Score Dream vs Healthy Controls
# ============================================================
cat("\n=== ANALYSIS 2: NAS Stages vs Healthy Controls ===\n")

NAS_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")

meta_nas_all <- meta[
  dataset %in% NAS_DATASETS &
  sample_id %in% pass_samples
]

# Control group: condition == "Control"
# Disease groups: samples with nas_score >= 1
# Exclude: NAS0 non-control samples (some patients score NAS=0 but are not healthy)
meta_ctrl_nas <- meta_nas_all[condition == "Control"]
meta_ctrl_nas[, grp := "Ctrl"]

meta_nas_dis <- meta_nas_all[
  !is.na(nas_score) &
  nas_score >= 1 &
  condition != "Control"
]
# Collapse NAS >= 7 into 7
meta_nas_dis[, nas_grp := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]
meta_nas_dis[, grp := paste0("NAS", nas_grp)]

meta_nas <- rbind(meta_ctrl_nas, meta_nas_dis, fill = TRUE)
cat("Reference (Control):", nrow(meta_ctrl_nas), "\n")
print(meta_ctrl_nas[, .N, by = dataset][order(dataset)])

nas_dist <- meta_nas[, .N, by = grp][order(grp)]
cat("\nSample distribution (group):\n")
print(nas_dist)
fwrite(nas_dist, file.path(SIGS, "nas_stage_vs_ctrl_sample_sizes.csv"))

# Subset DGE
idx_nas  <- colnames(dge) %in% meta_nas$sample_id
dge_nas  <- dge[, idx_nas]
cat("DGE subset:", ncol(dge_nas), "samples\n")

idx_order   <- match(colnames(dge_nas), meta_nas$sample_id)
matched_nas <- meta_nas[idx_order]
inferred_sex_nas <- meta_mat$inferred_sex[match(colnames(dge_nas), meta_mat$sample_id)]

nas_level_order <- c("Ctrl", paste0("NAS", 1:7))
info_nas <- data.frame(
  grp          = factor(matched_nas$grp, levels = nas_level_order),
  dataset      = factor(matched_nas$dataset),
  inferred_sex = factor(inferred_sex_nas),
  row.names    = colnames(dge_nas),
  stringsAsFactors = FALSE
)

cat("  grp level counts:\n")
print(table(info_nas$grp, useNA = "always"))

form_nas <- ~ grp + inferred_sex + (1 | dataset)
cat("\nDream formula:", deparse(form_nas), "\n")

nas_results <- run_stage_dream(
  dge_nas, info_nas,
  form   = form_nas,
  coef_prefix   = "grp",
  stage_levels  = paste0("NAS", 1:7),
  contrast_fmt  = "%s_vs_Ctrl",
  param  = param
)
setnames(nas_results, "adj.P.Val", "padj", skip_absent = TRUE)
fwrite(nas_results, file.path(SIGS, "nas_stage_vs_ctrl_dream.csv"))
cat("\nNAS vs control dream results saved:", nrow(nas_results), "rows\n")

cat("\n=== 14b_stage_vs_healthy_dream completed:", as.character(Sys.time()), "===\n")
