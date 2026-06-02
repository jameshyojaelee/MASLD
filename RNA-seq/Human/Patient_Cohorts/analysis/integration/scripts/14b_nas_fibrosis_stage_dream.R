#!/usr/bin/env Rscript
# 14b_nas_fibrosis_stage_dream.R
# Per-NAS-score and per-fibrosis-stage dream mega-analyses for Figure 2
# Produces: nas_score_dream.csv, fibrosis_stage_dream.csv
#
# Usage: Rscript 14b_nas_fibrosis_stage_dream.R
# SLURM: bigmem, 16 CPU, 200GB RAM, ~2-3h

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
# (Required to prevent findbars-related errors on this R/variancePartition version)
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

cat("=== 14b: Per-NAS-Score & Per-Fibrosis-Stage Dream ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Setup parallel backend ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load merged counts and metadata ---
cat("Loading merged counts...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))

cat("Loading matched metadata...\n")
meta_matched <- readRDS(file.path(RDIR, "integration/meta_matched.rds"))

# Load unified metadata for NAS/fibrosis columns
meta_unified <- fread(file.path(INT, "metadata/unified_metadata.csv"))
cat("Unified metadata:", nrow(meta_unified), "samples,", ncol(meta_unified), "columns\n")

# Load QC report for pass_technical filter
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]
cat("QC-passing samples (pass_technical):", length(pass_samples), "\n\n")

# ============================================================
# ANALYSIS 1: Per-NAS-Score Dream Mega-Analysis
# ============================================================
cat("=== ANALYSIS 1: Per-NAS-Score Dream ===\n")

# Datasets with per-sample NAS scores
NAS_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")

# Filter to NAS-annotated samples that pass QC
meta_nas <- meta_unified[
  dataset %in% NAS_DATASETS &
  !is.na(nas_score) &
  sample_id %in% pass_samples
]
cat("NAS-annotated QC-passing samples:", nrow(meta_nas), "\n")

# Print sample distribution per NAS score
nas_dist <- meta_nas[, .N, by = nas_score][order(nas_score)]
cat("\nSample distribution by NAS score:\n")
print(nas_dist)

# Collapse NAS >= 7 into "7+" if any level has n < 10
meta_nas[, nas_group := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]
nas_group_dist <- meta_nas[, .N, by = nas_group][order(nas_group)]
cat("\nAfter collapsing 7+:\n")
print(nas_group_dist)

# Warn if NAS=0 has very few samples
n_nas0 <- sum(meta_nas$nas_group == 0)
if (n_nas0 < 10) {
  cat("WARNING: NAS=0 has only", n_nas0, "samples — consider collapsing 0-1\n")
}

# Subset DGE to NAS samples
nas_samples <- meta_nas$sample_id
idx_nas <- colnames(dge) %in% nas_samples
dge_nas <- dge[, idx_nas]
cat("DGE subset:", ncol(dge_nas), "samples,", nrow(dge_nas), "genes\n")

# Verify intersection
cat("Sample overlap check:", sum(idx_nas), "of", length(nas_samples), "NAS samples found in DGE\n")

# Build info dataframe
matched_sex_nas <- meta_matched$inferred_sex[match(colnames(dge_nas), meta_matched$sample_id)]
matched_meta_nas <- meta_nas[match(colnames(dge_nas), meta_nas$sample_id)]

info_nas <- data.frame(
  nas_group   = factor(matched_meta_nas$nas_group),
  dataset     = factor(matched_meta_nas$dataset),
  inferred_sex = factor(matched_sex_nas),
  row.names   = colnames(dge_nas),
  stringsAsFactors = FALSE
)

# Set NAS=0 as reference level
info_nas$nas_group <- relevel(info_nas$nas_group, ref = "0")

# Filter by expression
keep_nas <- filterByExpr(dge_nas, group = info_nas$nas_group)
dge_nas <- dge_nas[keep_nas, , keep.lib.sizes = FALSE]
dge_nas <- calcNormFactors(dge_nas, method = "TMM")
cat("After expression filter:", nrow(dge_nas), "genes\n")

# Dream formula: NAS as factor + sex + random intercept per dataset
# Note: age excluded — GSE135251 (largest NAS dataset) lacks age data
form_nas <- ~ nas_group + inferred_sex + (1 | dataset)
cat("\nDream formula:", deparse(form_nas), "\n")

cat("Running voomWithDreamWeights...\n")
v_nas <- voomWithDreamWeights(dge_nas, form_nas, info_nas, BPPARAM = param)

cat("Running dream...\n")
fit_nas <- dream(v_nas, form_nas, info_nas, BPPARAM = param)
# DO NOT call eBayes() — dream() already computes moderated t-stats

# Extract results for each NAS level vs NAS=0
nas_levels <- levels(info_nas$nas_group)
nas_levels <- nas_levels[nas_levels != "0"]  # Skip reference

all_nas_results <- list()
for (lvl in nas_levels) {
  coef_name <- paste0("nas_group", lvl)
  if (coef_name %in% colnames(coef(fit_nas))) {
    tt <- topTable(fit_nas, coef = coef_name, number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)
    tt$nas_level <- as.integer(lvl)
    tt$contrast <- paste0("NAS", lvl, "_vs_NAS0")
    all_nas_results[[lvl]] <- as.data.table(tt)
    n_sig <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
    cat("  NAS", lvl, "vs NAS 0:", n_sig, "DEGs (padj<0.05, |LFC|>0.5)\n")
  } else {
    cat("  WARNING: Coefficient", coef_name, "not found in model\n")
  }
}

nas_results <- rbindlist(all_nas_results, fill = TRUE)
setnames(nas_results, "adj.P.Val", "padj", skip_absent = TRUE)
fwrite(nas_results, file.path(SIGS, "nas_score_dream.csv"))
cat("\nNAS dream results saved:", nrow(nas_results), "rows\n")

# Save sample size summary for figure annotation
fwrite(nas_group_dist, file.path(SIGS, "nas_score_sample_sizes.csv"))

# ============================================================
# ANALYSIS 2: Per-Fibrosis-Stage Dream Mega-Analysis
# ============================================================
cat("\n=== ANALYSIS 2: Per-Fibrosis-Stage Dream ===\n")

# Datasets with individual fibrosis staging (F0-F4)
# GSE213621 is excluded for clean per-stage contrasts: it has grouped
# fibrosis (F0F1, F2, F3F4) not individual stages.
# (PRJNA512027 was permanently removed from the pipeline 2026-05-15.)
# This leaves 6 datasets (~800-900 QC-passing samples).
FIB_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478",
                   "GSE193066", "GSE240729")

meta_fib <- meta_unified[
  dataset %in% FIB_DATASETS &
  !is.na(fibrosis_stage) &
  sample_id %in% pass_samples
]
cat("Fibrosis-annotated QC-passing samples:", nrow(meta_fib), "\n")

# Ensure fibrosis_stage is integer 0-4
meta_fib[, fib_stage := as.integer(fibrosis_stage)]
meta_fib <- meta_fib[fib_stage %in% 0:4]

# Print sample distribution
fib_dist <- meta_fib[, .N, by = fib_stage][order(fib_stage)]
cat("\nSample distribution by fibrosis stage:\n")
print(fib_dist)

# Subset DGE
fib_samples <- meta_fib$sample_id
idx_fib <- colnames(dge) %in% fib_samples
dge_fib <- dge[, idx_fib]
cat("DGE subset:", ncol(dge_fib), "samples,", nrow(dge_fib), "genes\n")

# Verify intersection
cat("Sample overlap check:", sum(idx_fib), "of", length(fib_samples), "fibrosis samples found in DGE\n")

# Build info dataframe
matched_sex_fib <- meta_matched$inferred_sex[match(colnames(dge_fib), meta_matched$sample_id)]
matched_meta_fib <- meta_fib[match(colnames(dge_fib), meta_fib$sample_id)]

info_fib <- data.frame(
  fib_stage    = factor(matched_meta_fib$fib_stage),
  dataset      = factor(matched_meta_fib$dataset),
  inferred_sex = factor(matched_sex_fib),
  row.names    = colnames(dge_fib),
  stringsAsFactors = FALSE
)

# Set F0 as reference
info_fib$fib_stage <- relevel(info_fib$fib_stage, ref = "0")

# Filter by expression
keep_fib <- filterByExpr(dge_fib, group = info_fib$fib_stage)
dge_fib <- dge_fib[keep_fib, , keep.lib.sizes = FALSE]
dge_fib <- calcNormFactors(dge_fib, method = "TMM")
cat("After expression filter:", nrow(dge_fib), "genes\n")

# Dream formula
form_fib <- ~ fib_stage + inferred_sex + (1 | dataset)
cat("\nDream formula:", deparse(form_fib), "\n")

cat("Running voomWithDreamWeights...\n")
v_fib <- voomWithDreamWeights(dge_fib, form_fib, info_fib, BPPARAM = param)

cat("Running dream...\n")
fit_fib <- dream(v_fib, form_fib, info_fib, BPPARAM = param)

# Extract results for each fibrosis stage vs F0
fib_levels <- levels(info_fib$fib_stage)
fib_levels <- fib_levels[fib_levels != "0"]

all_fib_results <- list()
for (lvl in fib_levels) {
  coef_name <- paste0("fib_stage", lvl)
  if (coef_name %in% colnames(coef(fit_fib))) {
    tt <- topTable(fit_fib, coef = coef_name, number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)
    tt$fib_stage <- as.integer(lvl)
    tt$contrast <- paste0("F", lvl, "_vs_F0")
    all_fib_results[[lvl]] <- as.data.table(tt)
    n_sig <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
    cat("  F", lvl, "vs F0:", n_sig, "DEGs (padj<0.05, |LFC|>0.5)\n")
  } else {
    cat("  WARNING: Coefficient", coef_name, "not found in model\n")
  }
}

fib_results <- rbindlist(all_fib_results, fill = TRUE)
setnames(fib_results, "adj.P.Val", "padj", skip_absent = TRUE)
fwrite(fib_results, file.path(SIGS, "fibrosis_stage_dream.csv"))
cat("\nFibrosis stage dream results saved:", nrow(fib_results), "rows\n")

# Save sample size summary
fwrite(fib_dist, file.path(SIGS, "fibrosis_stage_sample_sizes.csv"))

cat("\n=== 14b completed:", as.character(Sys.time()), "===\n")
