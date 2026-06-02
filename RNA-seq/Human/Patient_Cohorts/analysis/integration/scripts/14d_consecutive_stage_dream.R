#!/usr/bin/env Rscript
# 14d_consecutive_stage_dream.R
# Consecutive pairwise dream contrasts: NAS(k) vs NAS(k-1) and F(k) vs F(k-1)
# Uses treatment-coded dream + contrasts.fit() for consecutive differences
# Also computes progression landscape data (NAS × Fibrosis signature heatmap)
#
# Outputs:
#   nas_consecutive_dream.csv      - consecutive NAS contrast results
#   fibrosis_consecutive_dream.csv - consecutive fibrosis contrast results
#   progression_landscape_cells.csv  - per-cell mean signature scores
#   progression_landscape_samples.csv - per-sample signature scores
#   progression_sample_distribution.csv - NAS × Fibrosis sample counts
#   nas_signature_genes.csv, fibrosis_signature_genes.csv
#
# Usage: Rscript 14d_consecutive_stage_dream.R
# SLURM: bigmem, 16 CPU, 260GB RAM, ~3-4h

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(limma)
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

cat("=== 14d: Consecutive Stage Dream Contrasts ===\n")
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

meta_unified <- fread(file.path(INT, "metadata/unified_metadata.csv"))
cat("Unified metadata:", nrow(meta_unified), "samples,", ncol(meta_unified), "columns\n")

qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]
cat("QC-passing samples (pass_technical):", length(pass_samples), "\n\n")

# ============================================================
# ANALYSIS 1: NAS Consecutive Contrasts
# ============================================================
cat("=== ANALYSIS 1: NAS Consecutive Contrasts ===\n")

NAS_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")

meta_nas <- meta_unified[
  dataset %in% NAS_DATASETS &
  !is.na(nas_score) &
  sample_id %in% pass_samples
]
meta_nas[, nas_group := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]

nas_dist <- meta_nas[, .N, by = nas_group][order(nas_group)]
cat("NAS samples:", nrow(meta_nas), "\n")
cat("Distribution:\n")
print(nas_dist)

# Subset DGE
nas_samples <- meta_nas$sample_id
idx_nas <- colnames(dge) %in% nas_samples
dge_nas <- dge[, idx_nas]

matched_sex_nas <- meta_matched$inferred_sex[match(colnames(dge_nas), meta_matched$sample_id)]
matched_meta_nas <- meta_nas[match(colnames(dge_nas), meta_nas$sample_id)]

info_nas <- data.frame(
  nas_group    = factor(matched_meta_nas$nas_group),
  dataset      = factor(matched_meta_nas$dataset),
  inferred_sex = factor(matched_sex_nas),
  row.names    = colnames(dge_nas),
  stringsAsFactors = FALSE
)

# Treatment coding with NAS 0 as reference (default)
info_nas$nas_group <- relevel(info_nas$nas_group, ref = "0")
cat("\nNAS levels:", paste(levels(info_nas$nas_group), collapse = ", "), "\n")
cat("Using treatment coding (ref=0) + contrasts.fit for consecutive\n")

keep_nas <- filterByExpr(dge_nas, group = info_nas$nas_group)
dge_nas <- dge_nas[keep_nas, , keep.lib.sizes = FALSE]
dge_nas <- calcNormFactors(dge_nas, method = "TMM")
cat("After expression filter:", nrow(dge_nas), "genes\n")

form_nas <- ~ nas_group + inferred_sex + (1 | dataset)
cat("Dream formula:", deparse(form_nas), "\n")

cat("Running voomWithDreamWeights...\n")
v_nas <- voomWithDreamWeights(dge_nas, form_nas, info_nas, BPPARAM = param)

cat("Running dream...\n")
fit_nas <- dream(v_nas, form_nas, info_nas, BPPARAM = param)

cat("Model coefficients:", paste(colnames(coef(fit_nas)), collapse = ", "), "\n")

# Compute consecutive contrasts from treatment-coded dream results
# Dream coefficients: nas_group1 = NAS1-NAS0, nas_group2 = NAS2-NAS0, etc.
# First contrast (NAS1 vs NAS0) is already a coefficient — extract directly.
# Subsequent contrasts: NAS(k) vs NAS(k-1) computed from coefficient differences.
#
# For the first contrast we use topTable directly from the dream fit.
# For subsequent contrasts we compute logFC, t, and p from the coefficient
# covariance matrix (available per-gene via fit$cov.coefficients.list or
# the stdev.unscaled + sigma approach).
#
# Simplest robust approach: extract all cumulative-vs-NAS0 results from dream,
# then compute consecutive logFC differences and re-estimate significance using
# the per-gene residual variance.

nas_levels <- levels(info_nas$nas_group)
nas_levels_nz <- nas_levels[nas_levels != "0"]
coef_names <- colnames(coef(fit_nas))

# Extract all cumulative results first (these have proper dream moderated stats)
cumul_results <- list()
for (lvl in nas_levels_nz) {
  coef_name <- paste0("nas_group", lvl)
  tt <- topTable(fit_nas, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  cumul_results[[lvl]] <- as.data.table(tt)
}

# Build consecutive contrast results
# For NAS1 vs NAS0: use the nas_group1 coefficient directly (already consecutive)
# For NAS(k) vs NAS(k-1): logFC = logFC_k - logFC_(k-1)
#   SE = sqrt(SE_k^2 + SE_(k-1)^2 - 2*cov(k, k-1))
#   Since dream provides stdev.unscaled per coef, and cov.coefficients:
#   We use the simpler approach: compute from the raw coefficient matrix

beta_mat <- coef(fit_nas)  # genes × coefficients
su_mat   <- fit_nas$stdev.unscaled  # genes × coefficients (or single matrix)
sigma    <- fit_nas$sigma  # per-gene residual SD
df_resid <- fit_nas$df.residual  # per-gene residual df

# For dream, moderated df and sigma come from df.total and s2.post
if (!is.null(fit_nas$df.total)) {
  df_use <- fit_nas$df.total
} else {
  df_use <- df_resid
}
if (!is.null(fit_nas$s2.post)) {
  sigma_use <- sqrt(fit_nas$s2.post)
} else {
  sigma_use <- sigma
}

all_nas_consec <- list()
for (i in seq_along(nas_levels_nz)) {
  lvl <- nas_levels_nz[i]
  cname <- paste0("NAS", lvl, "_vs_NAS", ifelse(i == 1, "0", nas_levels_nz[i-1]))

  if (i == 1) {
    # First contrast = direct coefficient (NAS1 - NAS0)
    tt <- cumul_results[[lvl]]
    tt$contrast <- cname
    setnames(tt, "adj.P.Val", "padj", skip_absent = TRUE)
  } else {
    prev_lvl <- nas_levels_nz[i-1]
    curr_coef <- paste0("nas_group", lvl)
    prev_coef <- paste0("nas_group", prev_lvl)

    curr_idx <- which(coef_names == curr_coef)
    prev_idx <- which(coef_names == prev_coef)

    # logFC for consecutive contrast
    lfc <- beta_mat[, curr_idx] - beta_mat[, prev_idx]

    # SE: need covariance between coefficients
    # Dream stores cov.coefficients.list (per-gene) or cov.coefficients (shared)
    if (!is.null(fit_nas$cov.coefficients.list)) {
      # Per-gene covariance matrices
      se <- sapply(seq_len(nrow(beta_mat)), function(g) {
        V <- fit_nas$cov.coefficients.list[[g]]
        sqrt(V[curr_idx, curr_idx] + V[prev_idx, prev_idx] - 2*V[curr_idx, prev_idx])
      })
      se <- se * sigma_use
    } else if (!is.null(fit_nas$cov.coefficients)) {
      # Shared covariance structure
      V <- fit_nas$cov.coefficients
      se_unscaled <- sqrt(V[curr_idx, curr_idx] + V[prev_idx, prev_idx] - 2*V[curr_idx, prev_idx])
      se <- se_unscaled * sigma_use
    } else {
      # Fallback: assume independence (conservative)
      su_curr <- su_mat[, curr_idx]
      su_prev <- su_mat[, prev_idx]
      se <- sqrt(su_curr^2 + su_prev^2) * sigma_use
    }

    tstat <- lfc / se
    pval <- 2 * pt(abs(tstat), df = df_use, lower.tail = FALSE)
    padj <- p.adjust(pval, method = "BH")

    tt <- data.table(
      logFC = lfc,
      t = tstat,
      P.Value = pval,
      padj = padj,
      gene = rownames(beta_mat),
      contrast = cname
    )
  }

  all_nas_consec[[cname]] <- tt
  n_sig <- sum(tt$padj < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
  n_sig_loose <- sum(tt$padj < 0.1, na.rm = TRUE)
  mean_abs_lfc <- mean(abs(tt$logFC), na.rm = TRUE)
  cat("  ", cname, ": ", n_sig, " DEGs (padj<0.05, |LFC|>0.5)",
      " | ", n_sig_loose, " (padj<0.1)",
      " | mean|LFC|=", round(mean_abs_lfc, 3), "\n", sep = "")
}

nas_consec <- rbindlist(all_nas_consec, fill = TRUE)
setnames(nas_consec, "adj.P.Val", "padj", skip_absent = TRUE)
fwrite(nas_consec, file.path(SIGS, "nas_consecutive_dream.csv"))
cat("NAS consecutive results saved:", nrow(nas_consec), "rows\n\n")

# ============================================================
# ANALYSIS 2: Fibrosis Consecutive Contrasts
# ============================================================
cat("=== ANALYSIS 2: Fibrosis Consecutive Contrasts ===\n")

FIB_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478",
                   "GSE193066", "GSE240729")

meta_fib <- meta_unified[
  dataset %in% FIB_DATASETS &
  !is.na(fibrosis_stage) &
  sample_id %in% pass_samples
]
meta_fib[, fib_stage := as.integer(fibrosis_stage)]
meta_fib <- meta_fib[fib_stage %in% 0:4]

fib_dist <- meta_fib[, .N, by = fib_stage][order(fib_stage)]
cat("Fibrosis samples:", nrow(meta_fib), "\n")
cat("Distribution:\n")
print(fib_dist)

fib_samples <- meta_fib$sample_id
idx_fib <- colnames(dge) %in% fib_samples
dge_fib <- dge[, idx_fib]

matched_sex_fib <- meta_matched$inferred_sex[match(colnames(dge_fib), meta_matched$sample_id)]
matched_meta_fib <- meta_fib[match(colnames(dge_fib), meta_fib$sample_id)]

info_fib <- data.frame(
  fib_stage    = factor(matched_meta_fib$fib_stage),
  dataset      = factor(matched_meta_fib$dataset),
  inferred_sex = factor(matched_sex_fib),
  row.names    = colnames(dge_fib),
  stringsAsFactors = FALSE
)

info_fib$fib_stage <- relevel(info_fib$fib_stage, ref = "0")
cat("\nFibrosis levels:", paste(levels(info_fib$fib_stage), collapse = ", "), "\n")
cat("Using treatment coding (ref=0) + contrasts.fit for consecutive\n")

keep_fib <- filterByExpr(dge_fib, group = info_fib$fib_stage)
dge_fib <- dge_fib[keep_fib, , keep.lib.sizes = FALSE]
dge_fib <- calcNormFactors(dge_fib, method = "TMM")
cat("After expression filter:", nrow(dge_fib), "genes\n")

form_fib <- ~ fib_stage + inferred_sex + (1 | dataset)
cat("Dream formula:", deparse(form_fib), "\n")

cat("Running voomWithDreamWeights...\n")
v_fib <- voomWithDreamWeights(dge_fib, form_fib, info_fib, BPPARAM = param)

cat("Running dream...\n")
fit_fib <- dream(v_fib, form_fib, info_fib, BPPARAM = param)

cat("Model coefficients:", paste(colnames(coef(fit_fib)), collapse = ", "), "\n")

# Compute consecutive fibrosis contrasts (same approach as NAS above)
fib_levels <- levels(info_fib$fib_stage)
fib_levels_nz <- fib_levels[fib_levels != "0"]
coef_names_f <- colnames(coef(fit_fib))

# Extract cumulative results
cumul_results_f <- list()
for (lvl in fib_levels_nz) {
  coef_name <- paste0("fib_stage", lvl)
  tt <- topTable(fit_fib, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  cumul_results_f[[lvl]] <- as.data.table(tt)
}

beta_mat_f <- coef(fit_fib)
su_mat_f   <- fit_fib$stdev.unscaled
sigma_f    <- fit_fib$sigma

if (!is.null(fit_fib$df.total)) {
  df_use_f <- fit_fib$df.total
} else {
  df_use_f <- fit_fib$df.residual
}
if (!is.null(fit_fib$s2.post)) {
  sigma_use_f <- sqrt(fit_fib$s2.post)
} else {
  sigma_use_f <- sigma_f
}

all_fib_consec <- list()
for (i in seq_along(fib_levels_nz)) {
  lvl <- fib_levels_nz[i]
  cname <- paste0("F", lvl, "_vs_F", ifelse(i == 1, "0", fib_levels_nz[i-1]))

  if (i == 1) {
    tt <- cumul_results_f[[lvl]]
    tt$contrast <- cname
    setnames(tt, "adj.P.Val", "padj", skip_absent = TRUE)
  } else {
    prev_lvl <- fib_levels_nz[i-1]
    curr_coef <- paste0("fib_stage", lvl)
    prev_coef <- paste0("fib_stage", prev_lvl)

    curr_idx <- which(coef_names_f == curr_coef)
    prev_idx <- which(coef_names_f == prev_coef)

    lfc <- beta_mat_f[, curr_idx] - beta_mat_f[, prev_idx]

    if (!is.null(fit_fib$cov.coefficients.list)) {
      se <- sapply(seq_len(nrow(beta_mat_f)), function(g) {
        V <- fit_fib$cov.coefficients.list[[g]]
        sqrt(V[curr_idx, curr_idx] + V[prev_idx, prev_idx] - 2*V[curr_idx, prev_idx])
      })
      se <- se * sigma_use_f
    } else if (!is.null(fit_fib$cov.coefficients)) {
      V <- fit_fib$cov.coefficients
      se_unscaled <- sqrt(V[curr_idx, curr_idx] + V[prev_idx, prev_idx] - 2*V[curr_idx, prev_idx])
      se <- se_unscaled * sigma_use_f
    } else {
      su_curr <- su_mat_f[, curr_idx]
      su_prev <- su_mat_f[, prev_idx]
      se <- sqrt(su_curr^2 + su_prev^2) * sigma_use_f
    }

    tstat <- lfc / se
    pval <- 2 * pt(abs(tstat), df = df_use_f, lower.tail = FALSE)
    padj <- p.adjust(pval, method = "BH")

    tt <- data.table(
      logFC = lfc,
      t = tstat,
      P.Value = pval,
      padj = padj,
      gene = rownames(beta_mat_f),
      contrast = cname
    )
  }

  all_fib_consec[[cname]] <- tt
  n_sig <- sum(tt$padj < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
  n_sig_loose <- sum(tt$padj < 0.1, na.rm = TRUE)
  mean_abs_lfc <- mean(abs(tt$logFC), na.rm = TRUE)
  cat("  ", cname, ": ", n_sig, " DEGs (padj<0.05, |LFC|>0.5)",
      " | ", n_sig_loose, " (padj<0.1)",
      " | mean|LFC|=", round(mean_abs_lfc, 3), "\n", sep = "")
}

fib_consec <- rbindlist(all_fib_consec, fill = TRUE)
setnames(fib_consec, "adj.P.Val", "padj", skip_absent = TRUE)
fwrite(fib_consec, file.path(SIGS, "fibrosis_consecutive_dream.csv"))
cat("Fibrosis consecutive results saved:", nrow(fib_consec), "rows\n\n")

# ============================================================
# ANALYSIS 3: Progression Landscape Data
# ============================================================
cat("=== ANALYSIS 3: Progression Landscape Data ===\n")

# Samples with BOTH NAS and fibrosis scores (5 overlapping datasets)
BOTH_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")

meta_both <- meta_unified[
  dataset %in% BOTH_DATASETS &
  !is.na(nas_score) &
  !is.na(fibrosis_stage) &
  sample_id %in% pass_samples
]
meta_both[, nas_group := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]
meta_both[, fib_stage := as.integer(fibrosis_stage)]
meta_both <- meta_both[fib_stage %in% 0:4]
cat("Samples with both NAS + fibrosis:", nrow(meta_both), "\n")

# Sample distribution
sample_dist <- meta_both[, .N, by = .(nas_group, fib_stage)]
fwrite(sample_dist, file.path(SIGS, "progression_sample_distribution.csv"))
cat("Sample distribution:\n")
print(dcast(sample_dist, fib_stage ~ nas_group, value.var = "N", fill = 0))

# Define NAS and fibrosis signature genes from 14b cumulative results
nas_dream <- fread(file.path(SIGS, "nas_score_dream.csv"))
fib_dream <- fread(file.path(SIGS, "fibrosis_stage_dream.csv"))

# NAS signature: top 50 by min padj across all cumulative contrasts
nas_top <- nas_dream[, .(min_padj = min(padj, na.rm = TRUE),
                          max_abs_lfc = max(abs(logFC), na.rm = TRUE)), by = gene]
nas_top <- nas_top[order(min_padj)][1:50]
nas_sig_genes <- nas_top$gene
cat("NAS signature genes: top 50 by min padj\n")

# Fibrosis signature: top 50 by min padj across all cumulative contrasts
fib_top <- fib_dream[, .(min_padj = min(padj, na.rm = TRUE),
                           max_abs_lfc = max(abs(logFC), na.rm = TRUE)), by = gene]
fib_top <- fib_top[order(min_padj)][1:50]
fib_sig_genes <- fib_top$gene
cat("Fibrosis signature genes: top 50 by min padj\n")

# Save signature gene lists
fwrite(data.table(gene = nas_sig_genes, signature = "NAS"),
       file.path(SIGS, "nas_signature_genes.csv"))
fwrite(data.table(gene = fib_sig_genes, signature = "Fibrosis"),
       file.path(SIGS, "fibrosis_signature_genes.csv"))

# Compute CPM for landscape plot
both_samples <- meta_both$sample_id
idx_both <- colnames(dge) %in% both_samples
dge_both <- dge[, idx_both]
dge_both <- calcNormFactors(dge_both, method = "TMM")

cpm_mat <- cpm(dge_both, log = TRUE, prior.count = 1)  # log2(CPM + 1)
cat("CPM matrix:", nrow(cpm_mat), "genes x", ncol(cpm_mat), "samples\n")

# Strip Ensembl version for matching
rownames_base <- gsub("\\..*", "", rownames(cpm_mat))
nas_sig_base <- gsub("\\..*", "", nas_sig_genes)
fib_sig_base <- gsub("\\..*", "", fib_sig_genes)

nas_idx <- which(rownames_base %in% nas_sig_base)
fib_idx <- which(rownames_base %in% fib_sig_base)
cat("NAS signature genes found in CPM:", length(nas_idx), "\n")
cat("Fibrosis signature genes found in CPM:", length(fib_idx), "\n")

# Per-sample signature scores (mean expression of signature genes)
nas_sig_scores <- colMeans(cpm_mat[nas_idx, , drop = FALSE])
fib_sig_scores <- colMeans(cpm_mat[fib_idx, , drop = FALSE])

# Build per-sample data
landscape_samples <- data.table(
  sample_id = colnames(cpm_mat),
  nas_sig_score = nas_sig_scores,
  fib_sig_score = fib_sig_scores
)
matched_both <- meta_both[match(landscape_samples$sample_id, meta_both$sample_id)]
landscape_samples[, `:=`(
  nas_group = matched_both$nas_group,
  fib_stage = matched_both$fib_stage,
  dataset   = matched_both$dataset
)]

# Per-cell means for landscape heatmap
landscape_cells <- landscape_samples[, .(
  n_samples     = .N,
  nas_sig_mean  = mean(nas_sig_score),
  fib_sig_mean  = mean(fib_sig_score),
  nas_sig_sd    = sd(nas_sig_score),
  fib_sig_sd    = sd(fib_sig_score)
), by = .(nas_group, fib_stage)]

landscape_cells[, ratio := nas_sig_mean / fib_sig_mean]

fwrite(landscape_cells, file.path(SIGS, "progression_landscape_cells.csv"))
fwrite(landscape_samples, file.path(SIGS, "progression_landscape_samples.csv"))
cat("Landscape data saved:", nrow(landscape_cells), "cells\n")

# Signature correlation
cor_test <- cor.test(landscape_samples$nas_sig_score, landscape_samples$fib_sig_score)
cat("\nCorrelation between NAS & Fibrosis signatures:\n")
cat("  r =", round(cor_test$estimate, 3), "\n")
cat("  p =", format(cor_test$p.value, digits = 3), "\n")

cat("\n=== 14d completed:", as.character(Sys.time()), "===\n")
