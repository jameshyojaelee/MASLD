#!/usr/bin/env Rscript
# 25_deconv_attribution.R
# ---------------------------------------------------------------------------
# Strategy 2: Deconvolution Attribution Analysis
#
# Quantifies how much of the disease signal is driven by cell-type composition
# changes (Hepatocytes, Macrophages) versus intrinsic gene-expression changes.
#
# Steps:
#   1. Run unadjusted dream (no deconv covariates) as baseline
#   2. Compare adjusted (existing) vs unadjusted p-values → attribution scores
#   3. Classify genes: Hepatocyte_intrinsic, Composition_driven, Unmasked, NS
#   4. Run interaction DE: group_binary * Macrophages to find macrophage-
#      dependent disease genes
#   5. Generate 3 publication PDFs + 3 result CSVs
#
# Input:  merged_dge.rds, meta_matched.rds, MuSiC deconv proportions
# Output: results/causal_inference/deconv_attribution_scores.csv
#         results/causal_inference/interaction_de_results.csv
#         results/causal_inference/dream_results_unadjusted.csv
#         figures/deconv_venn.pdf
#         figures/deconv_attribution_heatmap.pdf
#         figures/deconv_interaction_volcano.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
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
  library(ggplot2)
})

# ============================================================
#  Paths
# ============================================================
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

OUTDIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference"
FIGDIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# T1.11 2026-04-22: Expand deconv covariate list from 2 → 16 cell types
# ---------------------------------------------------------------------------
# Team 1 §5 #5 audit: legacy adjusted model used only Hepatocytes + Macrophages
# (2 of 16 MuSiC cell types), biasing attribution to "hepatocyte-intrinsic /
# composition-driven" without accounting for T, B, NK, fibroblast, endothelial,
# etc. Expand to the full MuSiC output list so the adjustment is not partial.
#
# Team 1 §5 #2 also flagged threshold asymmetry: legacy used |LFC| > 0.5 for
# adjusted but |LFC| > 0.3 for unadjusted. As of the 0.3 -> 0.5 dream-LFC
# migration (LOO-CV stability; ~1.41x fold change), the canonical unadjusted
# threshold is 0.5, matching the adjusted model. Write symmetric-threshold
# sensitivity outputs alongside the legacy results (legacy kept for backward
# compatibility).
# ---------------------------------------------------------------------------
# Full 16-cell-type MuSiC column list (see MuSiC prop_weighted TSV header)
CT_COVAR_LIST_FULL <- c("Endothelial cells", "Hepatocytes", "Plasma cells", "T cells",
                        "Cholangiocytes", "Fibroblasts", "Macrophages",
                        "Circulating NK/NKT", "Resident NK",
                        "Mono+mono derived cells", "Basophils", "B cells",
                        "cDC1s", "cDC2s", "pDCs", "Neutrophils")
# Legacy (kept for compatibility when data missing): just Hep + Mac
CT_COVAR_LIST_LEGACY <- c("Hepatocytes", "Macrophages")
# Safe R variable names for the covariate columns in `info`
.t111_safe_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
CT_COVAR_VARS_FULL   <- .t111_safe_name(CT_COVAR_LIST_FULL)
CT_COVAR_VARS_LEGACY <- .t111_safe_name(CT_COVAR_LIST_LEGACY)

# ============================================================
#  Parallel setup
# ============================================================
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
if (ncpus > 1) {
  register(MulticoreParam(ncpus))
} else {
  register(SerialParam())
}

# ============================================================
#  1. Load data (same pattern as 05_dream_mega_analysis.R)
# ============================================================
cat("Loading merged DGE object...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Exclude datasets with irrecoverable confounds (same as script 05)
# GSE167523: no healthy controls.
# (PRJNA512027 — L0/S0 batch confounded with disease — was permanently removed
# from the pipeline 2026-05-15.)
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), "
")
cat("Samples for mega-analysis:", ncol(dge_mega), "\n")

# Reload metadata for inferred_sex
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# Load MuSiC deconvolution fractions
# T1.11: load ALL 16 cell types from MuSiC output, not just Hep+Mac.
cat("Loading MuSiC deconvolution fractions (all 16 cell types)...\n")
music_dir <- file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "Analysis/Deconvolution/results")
datasets <- unique(dge_mega$samples$dataset)
music_props <- lapply(datasets, function(ds) {
  prop_file <- file.path(music_dir, ds, paste0(ds, "_music_prop_weighted.tsv"))
  if (file.exists(prop_file)) {
    df <- read.table(prop_file, header = TRUE, sep = "\t", check.names = FALSE)
    df$sample_id <- rownames(df)
    # Keep any of CT_COVAR_LIST_FULL columns that are present
    present_cols <- intersect(CT_COVAR_LIST_FULL, colnames(df))
    keep <- c("sample_id", present_cols)
    return(df[, keep, drop = FALSE])
  } else {
    warning(paste("MuSiC proportions not found for dataset:", ds))
    return(NULL)
  }
})
music_dt <- rbindlist(music_props, fill = TRUE)

# Build a matched matrix of proportions, one column per cell type (full set)
# Use safe R variable names; impute missing cell types with dataset-wide median,
# and if entirely missing (cell type absent from MuSiC output), impute with 0.
matched_props <- matrix(NA_real_, nrow = ncol(dge_mega), ncol = length(CT_COVAR_LIST_FULL),
                        dimnames = list(colnames(dge_mega), CT_COVAR_VARS_FULL))
for (i in seq_along(CT_COVAR_LIST_FULL)) {
  ct_orig <- CT_COVAR_LIST_FULL[i]
  ct_safe <- CT_COVAR_VARS_FULL[i]
  if (ct_orig %in% colnames(music_dt)) {
    vec <- music_dt[[ct_orig]][match(colnames(dge_mega), music_dt$sample_id)]
    if (any(is.na(vec))) {
      med <- median(vec, na.rm = TRUE)
      if (is.na(med)) med <- 0
      vec[is.na(vec)] <- med
    }
    matched_props[, ct_safe] <- vec
  } else {
    # Cell type entirely absent in MuSiC outputs; fill with 0 (constant col
    # will be dropped by lm/dream automatically via alias detection).
    matched_props[, ct_safe] <- 0
  }
}

# Back-compat: matched_hepa / matched_macro still referenced by figures.
matched_hepa  <- matched_props[, .t111_safe_name("Hepatocytes")]
matched_macro <- matched_props[, .t111_safe_name("Macrophages")]

# Build design info: include all 16 covariates, plus dataset + sex + group
info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
for (ct_safe in CT_COVAR_VARS_FULL) {
  info[[ct_safe]] <- matched_props[, ct_safe]
}
# Back-compat aliases (used by plotting sections)
info$Hepatocytes <- matched_hepa
info$Macrophages <- matched_macro
rownames(info) <- colnames(dge_mega)

# Drop covariates that have zero variance (e.g. entirely absent from MuSiC)
# or median < 0.005 (rarely detected, adds noise + worsens conditioning).
# Additionally: MuSiC proportions are compositional (sum to 1), so including
# ALL cell types creates a rank-deficient design (last column = 1 - sum of
# others). Drop the largest cell type (Hepatocytes) as the compositional
# reference — other coefficients are then interpreted as changes relative to
# Hepatocyte fraction. This is the standard fix for compositional covariates.
.t111_varok <- sapply(CT_COVAR_VARS_FULL, function(v) var(info[[v]]) > 1e-6)
# Use mean > 0.001 threshold (MuSiC proportions are sparse; many immune cell
# types have median 0 but mean fraction ~0.01-0.1, carrying real signal).
.t111_median_ok <- sapply(CT_COVAR_VARS_FULL, function(v) mean(info[[v]]) > 0.001)
.t111_hep_safe <- .t111_safe_name("Hepatocytes")
CT_COVAR_VARS_KEEP <- CT_COVAR_VARS_FULL[.t111_varok & .t111_median_ok]
CT_COVAR_VARS_KEEP <- setdiff(CT_COVAR_VARS_KEEP, .t111_hep_safe)
# Defensive: ensure Macrophages is kept in the non-interaction set
# (it is needed to preserve the legacy interpretation + interaction later).
cat(sprintf("T1.11: %d of %d cell-type covariates included (var>1e-6, median>0.005, dropping %s as compositional reference).\n",
            length(CT_COVAR_VARS_KEEP), length(CT_COVAR_VARS_FULL), .t111_hep_safe))
cat("T1.11 covariates:", paste(CT_COVAR_VARS_KEEP, collapse = ", "), "\n")
# Scale the T1.11 covariates to unit SD to improve numerical conditioning in
# lme4. Write scaled versions into NEW columns (suffix _scaled) so the raw
# Hepatocytes / Macrophages columns remain for the legacy 2-CT fit + the
# interaction model (both expect raw proportions).
CT_COVAR_VARS_SCALED <- paste0(CT_COVAR_VARS_KEEP, "_s")
for (i in seq_along(CT_COVAR_VARS_KEEP)) {
  v <- CT_COVAR_VARS_KEEP[i]
  sdv <- sd(info[[v]])
  info[[ CT_COVAR_VARS_SCALED[i] ]] <- if (sdv > 0) {
    (info[[v]] - mean(info[[v]])) / sdv
  } else info[[v]]
}

cat("Group distribution:\n")
print(table(info$group_binary, info$dataset))

# ============================================================
#  2. Run ADJUSTED dream (with deconv covariates)
# NOTE: This model asks "after controlling for cell composition,
#       what disease-associated expression changes remain?"
#       Hepatocyte fraction falls in disease (composition-driven signal),
#       so DEGs that survive this adjustment are hepatocyte-intrinsic.
#       GSE213621 lacks MuSiC data → those samples get median-imputed fractions.
# ============================================================
cat("\n===== ADJUSTED DREAM (with ALL 16 deconv covariates, T1.11) =====\n")
# T1.11: full-CT formula using scaled covariates (exclude Hepatocytes as
# compositional reference — dropped upstream).
.t111_covar_terms <- paste(CT_COVAR_VARS_SCALED, collapse = " + ")
form_adj <- as.formula(sprintf(
  "~ group_binary + inferred_sex + %s + (1|dataset)", .t111_covar_terms
))
cat("Formula (T1.11 full-CT):", deparse(form_adj, width.cutoff = 500), "\n")

cat("Running voomWithDreamWeights (adjusted, full-CT)...\n")
v_adj <- suppressWarnings(voomWithDreamWeights(dge_mega, form_adj, info))

cat("Running dream() (adjusted, full-CT)...\n")
fit_adj <- suppressWarnings(dream(v_adj, form_adj, info))
# NOTE: do NOT call eBayes() after dream()

res_adj <- topTable(fit_adj, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res_adj$gene <- rownames(res_adj)
dt_adj <- as.data.table(res_adj)
setnames(dt_adj, "adj.P.Val", "padj")

cat("Adjusted (16ct) DEGs (padj < 0.05, |logFC| > 0.5):",
    nrow(dt_adj[padj < 0.05 & abs(logFC) > 0.5]), "\n")
cat("Adjusted (16ct) DEGs (padj < 0.05, |logFC| > 0.3, legacy sensitivity):",
    nrow(dt_adj[padj < 0.05 & abs(logFC) > 0.3]), "\n")

# Save adjusted results — NEW T1.11 output alongside the legacy file.
fwrite(dt_adj, file.path(OUTDIR, "dream_results_deconvolution_adjusted_16ct.csv"))
cat("Saved: dream_results_deconvolution_adjusted_16ct.csv (T1.11 full 16-CT fit)\n")

# ---- Legacy 2-covariate (Hep+Mac) adjusted fit, kept for backward compat ----
cat("\n[T1.11] Also running legacy 2-covariate (Hep+Mac) adjusted fit for reference...\n")
form_adj_legacy <- ~ group_binary + inferred_sex + Hepatocytes + Macrophages + (1|dataset)
v_adj_legacy <- suppressWarnings(voomWithDreamWeights(dge_mega, form_adj_legacy, info))
fit_adj_legacy <- suppressWarnings(dream(v_adj_legacy, form_adj_legacy, info))
res_adj_legacy <- topTable(fit_adj_legacy, coef = "group_binaryDisease",
                           number = Inf, sort.by = "none")
res_adj_legacy$gene <- rownames(res_adj_legacy)
dt_adj_legacy <- as.data.table(res_adj_legacy)
setnames(dt_adj_legacy, "adj.P.Val", "padj")
fwrite(dt_adj_legacy, file.path(OUTDIR, "dream_results_deconvolution_adjusted.csv"))
cat("Saved: dream_results_deconvolution_adjusted.csv (legacy Hep+Mac)\n")
# --- end T1.11 fit expansion ---

# ============================================================
#  3. Load UNADJUSTED dream results (primary, from script 05)
# ============================================================
cat("\nLoading unadjusted dream results (from 05_dream_mega_analysis.R)...\n")
# Use dream results with standard thresholds (padj < 0.05, |logFC| > 0.5)
ashr_path <- file.path(RDIR, "dream_results_ashr.csv")
if (file.exists(ashr_path)) {
  dt_unadj <- fread(ashr_path)
  # padj and logFC already exist as-is in the data
  dt_unadj[, dream_sig := !is.na(padj) & padj < 0.05 & abs(logFC) > 0.5]  # Migrated 0.3 -> 0.5 (LOO-CV stability)
  cat("  Using dream results (padj-based significance)\n")
  USE_ASHR_UNADJ <- TRUE
} else {
  dt_unadj <- fread(file.path(RDIR, "dream_results.csv"))
  USE_ASHR_UNADJ <- FALSE
  cat("  WARNING: dream_results_ashr.csv not found; falling back to dream_results.csv\n")
}

# Ensure both have the same gene set
common_genes <- intersect(dt_unadj$gene, dt_adj$gene)
cat("Common genes between adjusted/unadjusted:", length(common_genes), "\n")

dt_unadj_m <- dt_unadj[gene %in% common_genes]
dt_adj_m   <- dt_adj[gene %in% common_genes]

# Align by gene
setkey(dt_unadj_m, gene)
setkey(dt_adj_m, gene)

# ============================================================
#  4. Compute attribution scores and classify genes
# Terminology:
#   dt_adj_m   = ADJUSTED results (deconv-corrected, from this script)
#   dt_unadj_m = UNADJUSTED results (primary dream, from script 05)
#
#   Hepatocyte_intrinsic: sig in BOTH adjusted AND unadjusted
#     → gene survives deconv correction → genuine hepatocyte-driven signal
#   Composition_driven: sig in UNADJUSTED only (disappears when deconv added)
#     → gene is explained by cell-type composition change, not intrinsic expression
#   Unmasked: sig in ADJUSTED only (newly significant when deconv added)
#     → rare; gene was masked by composition variation, revealed by correction
# ============================================================
cat("\nComputing attribution scores...\n")

attrib <- data.table(
  gene         = dt_unadj_m$gene,
  logFC_adj    = dt_adj_m[dt_unadj_m$gene, logFC],
  logFC_unadj  = dt_unadj_m$logFC,
  padj_adj     = dt_adj_m[dt_unadj_m$gene, padj],
  padj_unadj   = dt_unadj_m$padj,
  AveExpr      = dt_unadj_m$AveExpr
)

# Attribution = 1 - (padj_adj / padj_unadj), clamped to [0, 1]
# High attribution → adjusted p-value much smaller than unadjusted → deconv revealed signal
# Negative attribution → adjusting made gene less significant → composition-driven
attrib[, attribution_raw := 1 - (padj_adj / padj_unadj)]
attrib[is.nan(attribution_raw), attribution_raw := 1]
attrib[attribution_raw > 1, attribution_raw := 1]
attrib[attribution_raw < 0, attribution_raw := 0]

# Gene classification: use padj for unadjusted (padj < 0.05, |logFC| > 0.5)
# Adjusted model (this script) still uses padj since it has no ashr step
attrib[, sig_adj   := padj_adj < 0.05 & abs(logFC_adj) > 0.5]
if (USE_ASHR_UNADJ) {
  # Unadjusted uses padj < 0.05, |logFC| > 0.5 (migrated 0.3 -> 0.5; LOO-CV stability)
  padj_vec <- dt_unadj_m[attrib$gene, padj]
  attrib[, sig_unadj := !is.na(padj_vec) & padj_vec < 0.05 & abs(logFC_unadj) > 0.5]
} else {
  attrib[, sig_unadj := padj_unadj < 0.05 & abs(logFC_unadj) > 0.5]
}

attrib[, category := fcase(
  sig_adj & sig_unadj,  "Hepatocyte_intrinsic",
  !sig_adj & sig_unadj, "Composition_driven",
  sig_adj & !sig_unadj, "Unmasked",
  default = "Not_significant"
)]

cat("\nGene classification summary:\n")
print(attrib[, .N, by = category][order(-N)])

# Save attribution scores (T1.11 16-CT-adjusted attribution + legacy 2-CT alias)
# The legacy filename (deconv_attribution_scores.csv) is still used as the
# primary consumer-facing output (read by Scripts 27a/36/38/80/202 + figures);
# per T1.11 it now reflects the full 16-CT adjustment rather than the
# previous Hep+Mac-only fit. The original 2-CT categorisation is still
# available via dream_results_deconvolution_adjusted.csv + a legacy-scored
# attribution table written below.
fwrite(attrib, file.path(OUTDIR, "deconv_attribution_scores.csv"))
cat("Saved: deconv_attribution_scores.csv (primary; 16-CT adjustment, T1.11)\n")
fwrite(attrib, file.path(OUTDIR, "deconv_attribution_scores_16ct.csv"))
cat("Saved: deconv_attribution_scores_16ct.csv (alias of primary, T1.11 tag)\n")

# ---- T1.11 symmetric-threshold sensitivity attribution ----
# Legacy call uses asymmetric thresholds (adj |LFC|>0.5, unadj |LFC|>0.3).
# Post-migration, canonical is symmetric 0.5/0.5; sensitivity also runs at 0.3/0.3
# for audit reproducibility (Team 1 §5 #2).
for (lfc_thr in c(0.3, 0.5)) {
  tmp <- copy(attrib)
  tmp[, sig_adj   := padj_adj < 0.05 & abs(logFC_adj)   > lfc_thr]
  tmp[, sig_unadj := padj_unadj < 0.05 & abs(logFC_unadj) > lfc_thr]
  tmp[, category := fcase(
    sig_adj & sig_unadj,  "Hepatocyte_intrinsic",
    !sig_adj & sig_unadj, "Composition_driven",
    sig_adj & !sig_unadj, "Unmasked",
    default = "Not_significant"
  )]
  fn <- sprintf("deconv_attribution_scores_16ct_symLFC_%.1f.csv", lfc_thr)
  fwrite(tmp, file.path(OUTDIR, fn))
  cat(sprintf("Saved: %s (symmetric |LFC|>%.1f; %d intrinsic, %d composition)\n",
              fn, lfc_thr,
              sum(tmp$category == "Hepatocyte_intrinsic"),
              sum(tmp$category == "Composition_driven")))
  rm(tmp)
}

# ---- T1.11 legacy 2-CT attribution (for backward-compat audit) ----
attrib_legacy <- data.table(
  gene        = dt_unadj_m$gene,
  logFC_adj   = dt_adj_legacy[dt_unadj_m$gene, logFC],
  logFC_unadj = dt_unadj_m$logFC,
  padj_adj    = dt_adj_legacy[dt_unadj_m$gene, padj],
  padj_unadj  = dt_unadj_m$padj,
  AveExpr     = dt_unadj_m$AveExpr
)
attrib_legacy[, attribution_raw := 1 - (padj_adj / padj_unadj)]
attrib_legacy[is.nan(attribution_raw), attribution_raw := 1]
attrib_legacy[attribution_raw > 1, attribution_raw := 1]
attrib_legacy[attribution_raw < 0, attribution_raw := 0]
attrib_legacy[, sig_adj   := padj_adj < 0.05 & abs(logFC_adj) > 0.5]
if (USE_ASHR_UNADJ) {
  padj_vec <- dt_unadj_m[attrib_legacy$gene, padj]
  attrib_legacy[, sig_unadj := !is.na(padj_vec) & padj_vec < 0.05 & abs(logFC_unadj) > 0.3]
} else {
  attrib_legacy[, sig_unadj := padj_unadj < 0.05 & abs(logFC_unadj) > 0.5]
}
attrib_legacy[, category := fcase(
  sig_adj & sig_unadj,  "Hepatocyte_intrinsic",
  !sig_adj & sig_unadj, "Composition_driven",
  sig_adj & !sig_unadj, "Unmasked",
  default = "Not_significant"
)]
fwrite(attrib_legacy, file.path(OUTDIR, "deconv_attribution_scores_legacy2ct.csv"))
cat("Saved: deconv_attribution_scores_legacy2ct.csv (legacy Hep+Mac only)\n")
# --- end T1.11 attribution outputs ---

# ============================================================
#  5. Interaction DE: group_binary * Macrophages
# ============================================================
cat("\n===== INTERACTION DE: Disease x Macrophages =====\n")
# T1.11: condition on all non-interaction cell-type covariates (scaled) plus
# raw Macrophages for the interaction. Drop Macrophages-scaled from the main
# effect to avoid collinearity with the raw Macrophages interaction term.
.t111_macro_safe   <- .t111_safe_name("Macrophages")
.t111_macro_scaled <- paste0(.t111_macro_safe, "_s")
.t111_inter_main   <- setdiff(CT_COVAR_VARS_SCALED, .t111_macro_scaled)
.t111_inter_terms  <- paste(.t111_inter_main, collapse = " + ")
form_inter <- as.formula(sprintf(
  "~ group_binary * Macrophages + inferred_sex + %s + (1|dataset)",
  .t111_inter_terms
))
cat("Formula (T1.11 full-CT interaction):", deparse(form_inter, width.cutoff = 500), "\n")

# NOTE: In MASLD, macrophage infiltration correlates with disease status
# (group_binary ~ Macrophages), which can cause rank deficiency in the
# fixed-effects matrix. Wrap in tryCatch to handle gracefully.
inter_result <- tryCatch({
  cat("Running voomWithDreamWeights (interaction)...\n")
  v_inter <- suppressWarnings(voomWithDreamWeights(dge_mega, form_inter, info))

  cat("Running dream() (interaction)...\n")
  fit_inter <- suppressWarnings(dream(v_inter, form_inter, info))

  # Extract the interaction coefficient
  inter_coef <- "group_binaryDisease:Macrophages"
  cat("Extracting coefficient:", inter_coef, "\n")
  res_inter <- topTable(fit_inter, coef = inter_coef, number = Inf, sort.by = "none")
  res_inter$gene <- rownames(res_inter)
  dt_inter <- as.data.table(res_inter)
  setnames(dt_inter, "adj.P.Val", "padj")
  cat("Interaction DEGs (padj < 0.05):", nrow(dt_inter[padj < 0.05]), "\n")
  cat("Interaction DEGs (padj < 0.05, |logFC| > 0.5):",
      nrow(dt_inter[padj < 0.05 & abs(logFC) > 0.5]), "\n")
  fwrite(dt_inter, file.path(OUTDIR, "interaction_de_results.csv"))
  cat("Saved: interaction_de_results.csv\n")
  dt_inter
}, error = function(e) {
  cat("WARNING: Interaction DE failed (likely rank-deficient model):", conditionMessage(e), "\n")
  cat("  Macrophage fraction is correlated with disease status in MASLD.\n")
  cat("  Skipping interaction DE — main deconvolution attribution results are unaffected.\n")
  NULL
})
dt_inter <- inter_result

# ============================================================
#  6. Load gene annotation for plot labels (GENCODE v49 cache)
# ============================================================
cat("\nLoading gene annotations for plot labels...\n")
CACHE_DIR <- file.path(INT, "results/gene_annotation")
gene_cache <- file.path(CACHE_DIR, "human_gencode_v49_genes.csv")

if (file.exists(gene_cache)) {
  gene_map <- fread(gene_cache)
  cat("  Loaded cached gene map:", nrow(gene_map), "genes\n")
} else {
  # Fallback: parse from GTF if cache not available
  HUMAN_GTF <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
  cat("  Parsing GTF for gene names...\n")
  awk_cmd <- paste0(
    "zcat ", HUMAN_GTF,
    " | awk -F'\\t' '$3==\"gene\"{",
    "  match($9,/gene_id \"([^\"]+)\"/,gid);",
    "  match($9,/gene_name \"([^\"]+)\"/,gn);",
    "  match($9,/gene_type \"([^\"]+)\"/,gt);",
    "  print gid[1],gn[1],gt[1]}' OFS='\\t'"
  )
  gene_map <- fread(cmd = awk_cmd, header = FALSE,
                    col.names = c("gene_id", "gene_name", "gene_type"))
  # Strip version suffix for matching
  gene_map[, gene_id_nover := sub("\\..*", "", gene_id)]
}

# Helper to get gene symbol from Ensembl ID
get_symbol <- function(ensembl_ids, map_dt) {
  # Try matching with version first, then without
  syms <- map_dt$gene_name[match(ensembl_ids, map_dt$gene_id)]
  missing <- is.na(syms)
  if (any(missing) && "gene_id_nover" %in% names(map_dt)) {
    ids_nover <- sub("\\..*", "", ensembl_ids[missing])
    syms[missing] <- map_dt$gene_name[match(ids_nover, map_dt$gene_id_nover)]
  }
  # Fallback to Ensembl ID if still missing
  syms[is.na(syms)] <- ensembl_ids[is.na(syms)]
  return(syms)
}

# ============================================================
#  7. FIGURE 1: Venn diagram — Adjusted vs Unadjusted DEGs
# ============================================================
cat("\n===== Generating Figure 1: Venn Diagram =====\n")

sig_adj_genes   <- attrib[sig_adj == TRUE, gene]
sig_unadj_genes <- attrib[sig_unadj == TRUE, gene]

n_adj_only   <- length(setdiff(sig_adj_genes, sig_unadj_genes))
n_unadj_only <- length(setdiff(sig_unadj_genes, sig_adj_genes))
n_both       <- length(intersect(sig_adj_genes, sig_unadj_genes))

cat("  Adjusted only (Unmasked):", n_adj_only, "\n")
cat("  Unadjusted only (Composition-driven):", n_unadj_only, "\n")
cat("  Both (Hepatocyte-intrinsic = survives deconv correction):", n_both, "\n")

# Euler-style Venn using ggplot (no VennDiagram dependency)
venn_df <- data.table(
  label = c(
    paste0("Adjusted only\n(Unmasked)\n", n_adj_only),
    paste0("Both\n(Intrinsic)\n", n_both),
    paste0("Unadjusted only\n(Composition)\n", n_unadj_only)
  ),
  x = c(-1.2, 0, 1.2),
  y = c(0, 0, 0),
  count = c(n_adj_only, n_both, n_unadj_only)
)

# Build a proper Venn via overlapping circles
theta <- seq(0, 2 * pi, length.out = 200)
circle_r <- 1.5
circle1 <- data.table(
  x = circle_r * cos(theta) - 0.7,
  y = circle_r * sin(theta),
  group = "Adjusted\n(w/ deconv)"
)
circle2 <- data.table(
  x = circle_r * cos(theta) + 0.7,
  y = circle_r * sin(theta),
  group = "Unadjusted\n(no deconv)"
)

p_venn <- ggplot() +
  geom_polygon(data = circle1, aes(x = x, y = y), fill = "#D81B60", alpha = 0.25, color = "#D81B60", linewidth = 0.8) +
  geom_polygon(data = circle2, aes(x = x, y = y), fill = "#1E88E5", alpha = 0.25, color = "#1E88E5", linewidth = 0.8) +
  annotate("text", x = -1.5, y = 0, label = n_adj_only, size = 6, fontface = "bold", color = "#D81B60") +
  annotate("text", x =  0.0, y = 0, label = n_both,     size = 6, fontface = "bold", color = "#333333") +
  annotate("text", x =  1.5, y = 0, label = n_unadj_only, size = 6, fontface = "bold", color = "#1E88E5") +
  annotate("text", x = -1.5, y = -0.5, label = "Unmasked", size = 3.5, color = "#D81B60") +
  annotate("text", x =  0.0, y = -0.5, label = "Intrinsic", size = 3.5, color = "#333333") +
  annotate("text", x =  1.5, y = -0.5, label = "Composition", size = 3.5, color = "#1E88E5") +
  annotate("text", x = -0.7, y = 1.8, label = "Adjusted\n(w/ deconv)", size = 3.5, fontface = "italic") +
  annotate("text", x =  0.7, y = 1.8, label = "Unadjusted\n(no deconv)", size = 3.5, fontface = "italic") +
  coord_fixed(xlim = c(-3, 3), ylim = c(-2, 2.5)) +
  labs(title = "DEG Overlap: Deconvolution-Adjusted vs Unadjusted",
       subtitle = ifelse(USE_ASHR_UNADJ,
                         "Unadjusted: padj < 0.05, |logFC| > 0.5; Adjusted: padj < 0.05, |logFC| > 0.5",
                         "Thresholds: padj < 0.05, |logFC| > 0.5")) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )

pdf(file.path(FIGDIR, "deconv_venn.pdf"), width = 7, height = 5)
print(p_venn)
dev.off()
cat("Saved: deconv_venn.pdf\n")

# ============================================================
#  8. FIGURE 2: Attribution Heatmap (top 50 genes)
# ============================================================
cat("\n===== Generating Figure 2: Attribution Heatmap =====\n")

# Select top 50 most significant genes (by min padj across adjusted/unadjusted)
attrib[, min_padj := pmin(padj_adj, padj_unadj, na.rm = TRUE)]
top50 <- head(attrib[order(min_padj)], 50)

# Annotate with gene symbols
top50[, gene_symbol := get_symbol(gene, gene_map)]

# Build heatmap data: logFC adjusted, logFC unadjusted, attribution
heat_dt <- melt(
  top50[, .(gene_symbol, logFC_adj, logFC_unadj, attribution_raw, category)],
  id.vars = c("gene_symbol", "category"),
  measure.vars = c("logFC_adj", "logFC_unadj", "attribution_raw"),
  variable.name = "metric",
  value.name = "value"
)

# Relabel metrics for clarity
heat_dt[, metric := fcase(
  metric == "logFC_adj",       "logFC (adjusted)",
  metric == "logFC_unadj",     "logFC (unadjusted)",
  metric == "attribution_raw", "Attribution score"
)]

# Order genes by adjusted logFC
gene_order <- top50[order(-logFC_adj), gene_symbol]
heat_dt[, gene_symbol := factor(gene_symbol, levels = rev(gene_order))]
heat_dt[, metric := factor(metric, levels = c("logFC (unadjusted)", "logFC (adjusted)", "Attribution score"))]

p_heat <- ggplot(heat_dt, aes(x = metric, y = gene_symbol, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#1E88E5", mid = "white", high = "#D81B60",
    midpoint = 0, name = "Value"
  ) +
  facet_grid(. ~ metric, scales = "free_x", space = "free_x") +
  labs(
    title = "Deconvolution Attribution: Top 50 DEGs",
    subtitle = "Comparing adjusted vs unadjusted disease effect",
    x = "", y = ""
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

pdf(file.path(FIGDIR, "deconv_attribution_heatmap.pdf"), width = 8, height = 10)
print(p_heat)
dev.off()
cat("Saved: deconv_attribution_heatmap.pdf\n")

# ============================================================
#  9. FIGURE 3: Interaction Volcano Plot
# ============================================================
cat("\n===== Generating Figure 3: Interaction Volcano =====\n")

if (!is.null(dt_inter)) {
dt_inter[, gene_symbol := get_symbol(gene, gene_map)]
dt_inter[, nlog10p := -log10(P.Value)]

# Significance categories for coloring
dt_inter[, sig_cat := fcase(
  padj < 0.05 & logFC >  0.5, "Up (Disease x Macrophage)",
  padj < 0.05 & logFC < -0.5, "Down (Disease x Macrophage)",
  default = "NS"
)]

# Label top genes
top_inter <- head(dt_inter[padj < 0.05][order(padj)], 20)

p_volc <- ggplot(dt_inter, aes(x = logFC, y = nlog10p, color = sig_cat)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_point(data = dt_inter[sig_cat != "NS"], alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c(
    "Up (Disease x Macrophage)"   = "#D81B60",
    "Down (Disease x Macrophage)" = "#1E88E5",
    "NS" = "grey70"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  labs(
    title = "Interaction DE: Disease x Macrophage Proportion",
    subtitle = paste0(
      "Genes with macrophage-dependent disease effects | ",
      nrow(dt_inter[padj < 0.05 & abs(logFC) > 0.5]), " significant"
    ),
    x = "Interaction logFC (Disease:Macrophages)",
    y = expression(-log[10](P)),
    color = ""
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Add labels for top genes if any are significant
if (nrow(top_inter) > 0) {
  p_volc <- p_volc +
    ggrepel::geom_text_repel(
      data = top_inter,
      aes(label = gene_symbol),
      size = 2.8,
      max.overlaps = 15,
      segment.color = "grey50",
      segment.size = 0.3,
      color = "black"
    )
}

pdf(file.path(FIGDIR, "deconv_interaction_volcano.pdf"), width = 8, height = 6)
print(p_volc)
dev.off()
cat("Saved: deconv_interaction_volcano.pdf\n")
} else {
  cat("Skipping interaction volcano plot (dt_inter is NULL from failed tryCatch).\n")
}

# ============================================================
#  Summary
# ============================================================
cat("\n===== DECONVOLUTION ATTRIBUTION ANALYSIS COMPLETE =====\n")
cat("Results saved to:", OUTDIR, "\n")
cat("Figures saved to:", FIGDIR, "\n")
cat("\nAttribution categories:\n")
print(attrib[, .N, by = category][order(-N)])
if (!is.null(dt_inter)) {
  cat("\nInteraction DE summary:\n")
  cat("  Total genes tested:", nrow(dt_inter), "\n")
  cat("  Significant interactions (padj < 0.05):", nrow(dt_inter[padj < 0.05]), "\n")
  cat("  Significant with |logFC| > 0.5:", nrow(dt_inter[padj < 0.05 & abs(logFC) > 0.5]), "\n")
} else {
  cat("\nInteraction DE: skipped (model failed).\n")
}
cat("\nDone:", format(Sys.time()), "\n")
