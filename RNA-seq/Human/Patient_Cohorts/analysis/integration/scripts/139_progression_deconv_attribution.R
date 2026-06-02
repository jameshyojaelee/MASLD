#!/usr/bin/env Rscript
# 139_progression_deconv_attribution.R
# ---------------------------------------------------------------------------
# Deconvolution Attribution for Progression Contrasts
#
# Re-runs the Script 25 attribution framework for progression contrasts
# instead of Disease-vs-Control:
#   C2: NASH vs NAFL (inflammation/steatohepatitis transition)
#   C3: Advanced (F3-F4) vs Early (F0-F2) fibrosis
#
# For each contrast, runs:
#   1. Unadjusted dream (contrast + sex + dataset random)
#   2. Adjusted dream (+ Hepatocytes + Macrophages fractions)
#   3. Attribution score: 1 - (padj_adj / padj_unadj), clamped [0,1]
#   4. Gene classification: Hepatocyte_intrinsic, Composition_driven,
#      Unmasked, Not_significant
#   5. Interaction model for macrophage-dependent effects
#   6. Comparison across C1/C2/C3 attribution profiles
#
# Input:
#   merged_counts_raw.rds, meta_matched.rds, MuSiC deconv proportions,
#   sample_qc_report.csv, modeling_metadata.csv
#
# Output (to results/progression/):
#   progression_deconv_attribution_c2.csv
#   progression_deconv_attribution_c3.csv
#   progression_deconv_adjusted_c2.csv
#   progression_deconv_adjusted_c3.csv
#   progression_deconv_interaction_c2.csv
#   progression_deconv_comparison.csv
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

# ============================================================
#  Paths
# ============================================================
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")
ODIR <- file.path(INT, "results/progression")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

# C1 attribution path (Script 25 output)
C1_ATTR_FILE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/deconv_attribution_scores.csv"

# MuSiC deconvolution directory
MUSIC_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/results"

# ============================================================
#  Parallel setup
# ============================================================
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
cat("Using", ncpus, "CPU cores\n")
BPPARAM <- if (ncpus > 1) MulticoreParam(ncpus, progressbar = TRUE) else SerialParam()

# ============================================================
#  1. Load data
# ============================================================
cat("Loading data...\n")
counts <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Load modeling metadata for NAS/fibrosis columns
modeling_meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"))

# Merge NAS/fibrosis annotations into meta
meta <- merge(meta, modeling_meta[, .(sample_id, fibrosis_stage, nas_score,
  diagnosis_harmonized, nas_group, fib_ge3, nas_ge5)],
  by = "sample_id", all.x = TRUE, suffixes = c("", ".mm"))

# Resolve duplicate columns -- prefer modeling_metadata values
for (col in c("fibrosis_stage", "nas_score", "diagnosis_harmonized")) {
  mm_col <- paste0(col, ".mm")
  if (mm_col %in% names(meta)) {
    na_mask <- is.na(meta[[col]]) | meta[[col]] == ""
    if (any(na_mask)) meta[[col]][na_mask] <- meta[[mm_col]][na_mask]
    meta[, (mm_col) := NULL]
  }
}

# Prepare sex covariate (inferred_sex preferred, fall back to annotated sex)
meta[, sex_covar := inferred_sex]
na_sex <- is.na(meta$sex_covar) | meta$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta)) {
  meta$sex_covar[na_sex] <- meta$sex[na_sex]
}
meta[, sex_covar := factor(sex_covar)]

cat(sprintf("Total QC-passing samples: %d\n", nrow(meta)))
cat(sprintf("Samples with fibrosis staging: %d\n", sum(!is.na(meta$fibrosis_stage))))
cat(sprintf("Samples with NAS score: %d\n", sum(!is.na(meta$nas_score))))
cat(sprintf("Samples with diagnosis_harmonized: %d\n",
  sum(!is.na(meta$diagnosis_harmonized) & meta$diagnosis_harmonized != "")))

# ============================================================
#  2. Load MuSiC deconvolution fractions
# ============================================================
cat("\nLoading MuSiC deconvolution fractions...\n")
all_datasets <- unique(meta$dataset)
music_props <- lapply(all_datasets, function(ds) {
  prop_file <- file.path(MUSIC_DIR, ds, paste0(ds, "_music_prop_weighted.tsv"))
  if (file.exists(prop_file)) {
    df <- read.table(prop_file, header = TRUE, sep = "\t", check.names = FALSE)
    df$sample_id <- rownames(df)
    return(df[, c("sample_id", "Hepatocytes", "Macrophages")])
  } else {
    warning(paste("MuSiC proportions not found for dataset:", ds))
    return(NULL)
  }
})
music_dt <- rbindlist(music_props, fill = TRUE)
cat(sprintf("  MuSiC fractions loaded for %d samples\n", nrow(music_dt)))

# ============================================================
#  Helper: attach MuSiC fractions to a metadata subset
# ============================================================
attach_music_fractions <- function(meta_sub) {
  matched_hepa  <- music_dt$Hepatocytes[match(meta_sub$sample_id, music_dt$sample_id)]
  matched_macro <- music_dt$Macrophages[match(meta_sub$sample_id, music_dt$sample_id)]

  # Impute missing fractions with median (same as Script 25)
  if (any(is.na(matched_hepa)))  matched_hepa[is.na(matched_hepa)]   <- median(matched_hepa, na.rm = TRUE)
  if (any(is.na(matched_macro))) matched_macro[is.na(matched_macro)] <- median(matched_macro, na.rm = TRUE)

  meta_sub[, Hepatocytes := matched_hepa]
  meta_sub[, Macrophages := matched_macro]

  n_imputed <- sum(is.na(music_dt$Hepatocytes[match(meta_sub$sample_id, music_dt$sample_id)]))
  cat(sprintf("  MuSiC: %d samples matched, %d imputed with median\n",
    nrow(meta_sub) - n_imputed, n_imputed))
  return(meta_sub)
}

# ============================================================
#  Helper: run dream for a contrast subset
# ============================================================
run_dream_for_contrast <- function(meta_sub, dge, form, info_cols) {
  info <- data.frame(
    row.names = meta_sub$sample_id,
    stringsAsFactors = FALSE
  )
  for (col in info_cols) {
    info[[col]] <- meta_sub[[col]]
  }

  # Factorize
  info$dataset <- factor(info$dataset)
  info$sex_covar <- factor(info$sex_covar)

  cat(sprintf("  Formula: %s\n", deparse(form)))
  cat(sprintf("  Running voomWithDreamWeights...\n"))
  v <- suppressWarnings(voomWithDreamWeights(dge, form, info, BPPARAM = BPPARAM))

  cat(sprintf("  Running dream()...\n"))
  fit <- suppressWarnings(dream(v, form, info, BPPARAM = BPPARAM))
  # NOTE: do NOT call eBayes() after dream()

  return(fit)
}

# ============================================================
#  Helper: compute attribution and classify genes
# ============================================================
compute_attribution <- function(dt_unadj, dt_adj, sig_threshold = 0.1) {
  # Ensure both have the same gene set
  common_genes <- intersect(dt_unadj$gene, dt_adj$gene)
  cat(sprintf("  Common genes: %d\n", length(common_genes)))

  dt_unadj_m <- dt_unadj[gene %in% common_genes]
  dt_adj_m   <- dt_adj[gene %in% common_genes]

  setkey(dt_unadj_m, gene)
  setkey(dt_adj_m, gene)

  attrib <- data.table(
    gene         = dt_unadj_m$gene,
    logFC_adj    = dt_adj_m[dt_unadj_m$gene, logFC],
    logFC_unadj  = dt_unadj_m$logFC,
    padj_adj     = dt_adj_m[dt_unadj_m$gene, padj],
    padj_unadj   = dt_unadj_m$padj,
    AveExpr      = dt_unadj_m$AveExpr
  )

  # Attribution = 1 - (padj_adj / padj_unadj), clamped to [0, 1]
  attrib[, attribution_raw := 1 - (padj_adj / padj_unadj)]
  attrib[attribution_raw > 1, attribution_raw := 1]
  attrib[attribution_raw < 0, attribution_raw := 0]

  # Gene classification (padj < sig_threshold)
  # NOTE: Using padj < 0.1 (not 0.05) to match the project-wide DEG threshold
  # from dream mega-analysis (see CLAUDE.md: "padj < 0.1; 16,333 genes").
  # Script 25 used padj < 0.05 & |logFC| > 0.5 — we use padj < 0.1 here for
  # consistency with the progression pipeline, which defines DEGs at padj < 0.1.
  attrib[, sig_adj   := padj_adj < sig_threshold]
  attrib[, sig_unadj := padj_unadj < sig_threshold]

  attrib[, category := fcase(
    sig_adj & sig_unadj,  "Hepatocyte_intrinsic",
    !sig_adj & sig_unadj, "Composition_driven",
    sig_adj & !sig_unadj, "Unmasked",
    default = "Not_significant"
  )]

  cat("\n  Gene classification summary:\n")
  print(attrib[, .N, by = category][order(-N)])

  return(attrib)
}

# ============================================================
#  Helper: filter datasets needing >= 2 samples per group
# ============================================================
filter_valid_datasets <- function(meta_sub, group_col, group_levels) {
  ds_counts <- meta_sub[, .(
    n1 = sum(.SD[[group_col]] == group_levels[1]),
    n2 = sum(.SD[[group_col]] == group_levels[2])
  ), by = dataset, .SDcols = group_col]
  valid_ds <- ds_counts[n1 >= 2 & n2 >= 2, dataset]
  cat(sprintf("  Valid datasets (>= 2 per group): %d (%s)\n",
    length(valid_ds), paste(valid_ds, collapse = ", ")))
  return(meta_sub[dataset %in% valid_ds])
}

# ============================================================
#  Helper: build DGEList from counts + metadata
# ============================================================
build_dge <- function(meta_sub, group_col) {
  keep_samples <- intersect(meta_sub$sample_id, colnames(counts))
  meta_sub <- meta_sub[sample_id %in% keep_samples]
  dge <- DGEList(counts = counts[, keep_samples])
  dge <- calcNormFactors(dge, method = "TMM")
  keep_genes <- filterByExpr(dge, group = meta_sub[[group_col]])
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  cat(sprintf("  DGE: %d samples, %d genes\n", ncol(dge), nrow(dge)))
  return(list(dge = dge, meta = meta_sub))
}

# ============================================================
#  C2: NASH vs NAFL Deconvolution Attribution
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  C2: NASH vs NAFL — DECONVOLUTION ATTRIBUTION\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Build C2 subset: NAFL + NASH (Borderline grouped with NASH, per Script 13)
meta_c2 <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline")]
meta_c2[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
meta_c2[, nafl_nash := factor(nafl_nash, levels = c("NAFL", "NASH"))]

cat("C2 sample breakdown:\n")
print(meta_c2[, .N, by = .(dataset, nafl_nash)][order(dataset, nafl_nash)])

meta_c2 <- filter_valid_datasets(meta_c2, "nafl_nash", c("NAFL", "NASH"))
meta_c2 <- attach_music_fractions(meta_c2)

# Build DGEList
c2_built <- build_dge(meta_c2, "nafl_nash")
dge_c2   <- c2_built$dge
meta_c2  <- c2_built$meta

# --- C2 UNADJUSTED dream ---
cat("\n--- C2 UNADJUSTED (no deconv covariates) ---\n")
form_c2_unadj <- ~ nafl_nash + sex_covar + (1 | dataset)

info_c2_unadj <- data.frame(
  nafl_nash = factor(meta_c2$nafl_nash, levels = c("NAFL", "NASH")),
  sex_covar = factor(meta_c2$sex_covar),
  dataset   = factor(meta_c2$dataset),
  row.names = meta_c2$sample_id,
  stringsAsFactors = FALSE
)

cat(sprintf("  Formula: %s\n", deparse(form_c2_unadj)))
cat("  Running voomWithDreamWeights (unadjusted)...\n")
v_c2_unadj <- suppressWarnings(voomWithDreamWeights(dge_c2, form_c2_unadj, info_c2_unadj, BPPARAM = BPPARAM))
cat("  Running dream() (unadjusted)...\n")
fit_c2_unadj <- suppressWarnings(dream(v_c2_unadj, form_c2_unadj, info_c2_unadj, BPPARAM = BPPARAM))

res_c2_unadj <- topTable(fit_c2_unadj, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
res_c2_unadj$gene <- rownames(res_c2_unadj)
dt_c2_unadj <- as.data.table(res_c2_unadj)
setnames(dt_c2_unadj, "adj.P.Val", "padj")

cat(sprintf("  Unadjusted C2 DEGs (padj<0.1): %d\n", nrow(dt_c2_unadj[padj < 0.1])))
cat(sprintf("  Unadjusted C2 DEGs (padj<0.05): %d\n", nrow(dt_c2_unadj[padj < 0.05])))

# --- C2 ADJUSTED dream (with deconv covariates) ---
cat("\n--- C2 ADJUSTED (with Hepatocytes + Macrophages) ---\n")
form_c2_adj <- ~ nafl_nash + sex_covar + Hepatocytes + Macrophages + (1 | dataset)

info_c2_adj <- data.frame(
  nafl_nash   = factor(meta_c2$nafl_nash, levels = c("NAFL", "NASH")),
  sex_covar   = factor(meta_c2$sex_covar),
  Hepatocytes = meta_c2$Hepatocytes,
  Macrophages = meta_c2$Macrophages,
  dataset     = factor(meta_c2$dataset),
  row.names   = meta_c2$sample_id,
  stringsAsFactors = FALSE
)

cat(sprintf("  Formula: %s\n", deparse(form_c2_adj)))
cat("  Running voomWithDreamWeights (adjusted)...\n")
v_c2_adj <- suppressWarnings(voomWithDreamWeights(dge_c2, form_c2_adj, info_c2_adj, BPPARAM = BPPARAM))
cat("  Running dream() (adjusted)...\n")
fit_c2_adj <- suppressWarnings(dream(v_c2_adj, form_c2_adj, info_c2_adj, BPPARAM = BPPARAM))

res_c2_adj <- topTable(fit_c2_adj, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
res_c2_adj$gene <- rownames(res_c2_adj)
dt_c2_adj <- as.data.table(res_c2_adj)
setnames(dt_c2_adj, "adj.P.Val", "padj")

cat(sprintf("  Adjusted C2 DEGs (padj<0.1): %d\n", nrow(dt_c2_adj[padj < 0.1])))
cat(sprintf("  Adjusted C2 DEGs (padj<0.05): %d\n", nrow(dt_c2_adj[padj < 0.05])))

# Save adjusted results
fwrite(dt_c2_adj, file.path(ODIR, "progression_deconv_adjusted_c2.csv"))
cat("  Saved: progression_deconv_adjusted_c2.csv\n")

# --- C2 Attribution ---
cat("\n--- C2 ATTRIBUTION ---\n")
attrib_c2 <- compute_attribution(dt_c2_unadj, dt_c2_adj, sig_threshold = 0.1)
attrib_c2[, contrast := "C2_NASH_vs_NAFL"]
fwrite(attrib_c2, file.path(ODIR, "progression_deconv_attribution_c2.csv"))
cat("  Saved: progression_deconv_attribution_c2.csv\n")

# --- C2 Interaction: nafl_nash * Macrophages ---
cat("\n--- C2 INTERACTION: NAFL_NASH x Macrophages ---\n")
form_c2_inter <- ~ nafl_nash * Macrophages + sex_covar + Hepatocytes + (1 | dataset)

inter_c2_result <- tryCatch({
  cat(sprintf("  Formula: %s\n", deparse(form_c2_inter)))

  info_c2_inter <- data.frame(
    nafl_nash   = factor(meta_c2$nafl_nash, levels = c("NAFL", "NASH")),
    sex_covar   = factor(meta_c2$sex_covar),
    Hepatocytes = meta_c2$Hepatocytes,
    Macrophages = meta_c2$Macrophages,
    dataset     = factor(meta_c2$dataset),
    row.names   = meta_c2$sample_id,
    stringsAsFactors = FALSE
  )

  cat("  Running voomWithDreamWeights (interaction)...\n")
  v_c2_inter <- suppressWarnings(voomWithDreamWeights(dge_c2, form_c2_inter, info_c2_inter, BPPARAM = BPPARAM))
  cat("  Running dream() (interaction)...\n")
  fit_c2_inter <- suppressWarnings(dream(v_c2_inter, form_c2_inter, info_c2_inter, BPPARAM = BPPARAM))

  inter_coef <- "nafl_nashNASH:Macrophages"
  cat(sprintf("  Extracting coefficient: %s\n", inter_coef))
  res_c2_inter <- topTable(fit_c2_inter, coef = inter_coef, number = Inf, sort.by = "none")
  res_c2_inter$gene <- rownames(res_c2_inter)
  dt_c2_inter <- as.data.table(res_c2_inter)
  setnames(dt_c2_inter, "adj.P.Val", "padj")
  dt_c2_inter[, contrast := "C2_NASH_vs_NAFL"]

  cat(sprintf("  Interaction DEGs (padj<0.1): %d\n", nrow(dt_c2_inter[padj < 0.1])))
  cat(sprintf("  Interaction DEGs (padj<0.05): %d\n", nrow(dt_c2_inter[padj < 0.05])))

  fwrite(dt_c2_inter, file.path(ODIR, "progression_deconv_interaction_c2.csv"))
  cat("  Saved: progression_deconv_interaction_c2.csv\n")
  dt_c2_inter
}, error = function(e) {
  cat(sprintf("  WARNING: Interaction DE failed (likely rank-deficient): %s\n", conditionMessage(e)))
  cat("  Macrophage proportion correlates with NASH status — model may be rank-deficient.\n")
  cat("  Skipping interaction DE.\n")
  NULL
})

gc()

# ============================================================
#  C3: Advanced vs Early Fibrosis Deconvolution Attribution
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  C3: ADVANCED (F3-F4) vs EARLY (F0-F2) FIBROSIS — DECONVOLUTION ATTRIBUTION\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Build C3 subset: F0-F2 vs F3-F4 (consistent with Script 130)
meta_c3 <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c3[, fib_group := fifelse(fibrosis_stage >= 3, "Advanced", "Early")]
meta_c3[, fib_group := factor(fib_group, levels = c("Early", "Advanced"))]

cat("C3 sample breakdown:\n")
print(meta_c3[, .N, by = .(dataset, fib_group)][order(dataset, fib_group)])

meta_c3 <- filter_valid_datasets(meta_c3, "fib_group", c("Early", "Advanced"))
meta_c3 <- attach_music_fractions(meta_c3)

# Build DGEList
c3_built <- build_dge(meta_c3, "fib_group")
dge_c3   <- c3_built$dge
meta_c3  <- c3_built$meta

# --- C3 UNADJUSTED dream ---
cat("\n--- C3 UNADJUSTED (no deconv covariates) ---\n")
form_c3_unadj <- ~ fib_group + sex_covar + (1 | dataset)

info_c3_unadj <- data.frame(
  fib_group = factor(meta_c3$fib_group, levels = c("Early", "Advanced")),
  sex_covar = factor(meta_c3$sex_covar),
  dataset   = factor(meta_c3$dataset),
  row.names = meta_c3$sample_id,
  stringsAsFactors = FALSE
)

cat(sprintf("  Formula: %s\n", deparse(form_c3_unadj)))
cat("  Running voomWithDreamWeights (unadjusted)...\n")
v_c3_unadj <- suppressWarnings(voomWithDreamWeights(dge_c3, form_c3_unadj, info_c3_unadj, BPPARAM = BPPARAM))
cat("  Running dream() (unadjusted)...\n")
fit_c3_unadj <- suppressWarnings(dream(v_c3_unadj, form_c3_unadj, info_c3_unadj, BPPARAM = BPPARAM))

res_c3_unadj <- topTable(fit_c3_unadj, coef = "fib_groupAdvanced", number = Inf, sort.by = "none")
res_c3_unadj$gene <- rownames(res_c3_unadj)
dt_c3_unadj <- as.data.table(res_c3_unadj)
setnames(dt_c3_unadj, "adj.P.Val", "padj")

cat(sprintf("  Unadjusted C3 DEGs (padj<0.1): %d\n", nrow(dt_c3_unadj[padj < 0.1])))
cat(sprintf("  Unadjusted C3 DEGs (padj<0.05): %d\n", nrow(dt_c3_unadj[padj < 0.05])))

# --- C3 ADJUSTED dream (with deconv covariates) ---
cat("\n--- C3 ADJUSTED (with Hepatocytes + Macrophages) ---\n")
form_c3_adj <- ~ fib_group + sex_covar + Hepatocytes + Macrophages + (1 | dataset)

info_c3_adj <- data.frame(
  fib_group   = factor(meta_c3$fib_group, levels = c("Early", "Advanced")),
  sex_covar   = factor(meta_c3$sex_covar),
  Hepatocytes = meta_c3$Hepatocytes,
  Macrophages = meta_c3$Macrophages,
  dataset     = factor(meta_c3$dataset),
  row.names   = meta_c3$sample_id,
  stringsAsFactors = FALSE
)

cat(sprintf("  Formula: %s\n", deparse(form_c3_adj)))
cat("  Running voomWithDreamWeights (adjusted)...\n")
v_c3_adj <- suppressWarnings(voomWithDreamWeights(dge_c3, form_c3_adj, info_c3_adj, BPPARAM = BPPARAM))
cat("  Running dream() (adjusted)...\n")
fit_c3_adj <- suppressWarnings(dream(v_c3_adj, form_c3_adj, info_c3_adj, BPPARAM = BPPARAM))

res_c3_adj <- topTable(fit_c3_adj, coef = "fib_groupAdvanced", number = Inf, sort.by = "none")
res_c3_adj$gene <- rownames(res_c3_adj)
dt_c3_adj <- as.data.table(res_c3_adj)
setnames(dt_c3_adj, "adj.P.Val", "padj")

cat(sprintf("  Adjusted C3 DEGs (padj<0.1): %d\n", nrow(dt_c3_adj[padj < 0.1])))
cat(sprintf("  Adjusted C3 DEGs (padj<0.05): %d\n", nrow(dt_c3_adj[padj < 0.05])))

# Save adjusted results
fwrite(dt_c3_adj, file.path(ODIR, "progression_deconv_adjusted_c3.csv"))
cat("  Saved: progression_deconv_adjusted_c3.csv\n")

# --- C3 Attribution ---
cat("\n--- C3 ATTRIBUTION ---\n")
attrib_c3 <- compute_attribution(dt_c3_unadj, dt_c3_adj, sig_threshold = 0.1)
attrib_c3[, contrast := "C3_Advanced_vs_Early"]
fwrite(attrib_c3, file.path(ODIR, "progression_deconv_attribution_c3.csv"))
cat("  Saved: progression_deconv_attribution_c3.csv\n")

gc()

# ============================================================
#  Comparison: C1 vs C2 vs C3 Attribution Profiles
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  COMPARISON: C1 vs C2 vs C3 ATTRIBUTION\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Load C1 (Disease vs Control, from Script 25)
if (file.exists(C1_ATTR_FILE)) {
  attrib_c1 <- fread(C1_ATTR_FILE)
  attrib_c1[, contrast := "C1_Disease_vs_Control"]
  cat(sprintf("  C1 attribution loaded: %d genes\n", nrow(attrib_c1)))
} else {
  cat("  WARNING: C1 attribution not found — comparison will only include C2 and C3.\n")
  attrib_c1 <- NULL
}

# Standardize column names across all contrasts
# C1 (Script 25) may have slightly different column layout; align
common_cols <- c("gene", "logFC_adj", "logFC_unadj", "padj_adj", "padj_unadj",
                 "AveExpr", "attribution_raw", "sig_adj", "sig_unadj", "category", "contrast")

format_attrib <- function(dt, expected_cols) {
  # Keep only the common columns that exist
  available <- intersect(expected_cols, names(dt))
  return(dt[, ..available])
}

attrib_c2_clean <- format_attrib(attrib_c2, common_cols)
attrib_c3_clean <- format_attrib(attrib_c3, common_cols)

if (!is.null(attrib_c1)) {
  attrib_c1_clean <- format_attrib(attrib_c1, common_cols)
  all_attrib <- rbindlist(list(attrib_c1_clean, attrib_c2_clean, attrib_c3_clean),
    use.names = TRUE, fill = TRUE)
} else {
  all_attrib <- rbindlist(list(attrib_c2_clean, attrib_c3_clean),
    use.names = TRUE, fill = TRUE)
}

# Summary table: per-contrast classification counts
comparison_summary <- all_attrib[, .(
  n_genes = .N,
  n_Hepatocyte_intrinsic = sum(category == "Hepatocyte_intrinsic"),
  n_Composition_driven   = sum(category == "Composition_driven"),
  n_Unmasked             = sum(category == "Unmasked"),
  n_Not_significant      = sum(category == "Not_significant"),
  pct_intrinsic          = round(100 * sum(category == "Hepatocyte_intrinsic") / .N, 1),
  pct_composition        = round(100 * sum(category == "Composition_driven") / .N, 1),
  pct_unmasked           = round(100 * sum(category == "Unmasked") / .N, 1),
  mean_attribution_sig   = round(mean(attribution_raw[sig_unadj == TRUE], na.rm = TRUE), 3),
  median_attribution_sig = round(median(attribution_raw[sig_unadj == TRUE], na.rm = TRUE), 3)
), by = contrast]

cat("\n  Attribution comparison summary:\n")
print(comparison_summary)

# Intrinsic fraction among DEGs only
intrinsic_among_degs <- all_attrib[sig_unadj == TRUE, .(
  n_degs = .N,
  n_intrinsic = sum(category == "Hepatocyte_intrinsic"),
  pct_intrinsic = round(100 * sum(category == "Hepatocyte_intrinsic") / .N, 1),
  n_composition = sum(category == "Composition_driven"),
  pct_composition = round(100 * sum(category == "Composition_driven") / .N, 1)
), by = contrast]

cat("\n  Among DEGs (sig_unadj == TRUE):\n")
print(intrinsic_among_degs)

# Per-gene comparison across contrasts (wide format)
# Find genes present in multiple contrasts
gene_contrast_wide <- dcast(
  all_attrib[, .(gene, contrast, category, attribution_raw)],
  gene ~ contrast,
  value.var = list("category", "attribution_raw"),
  fill = NA
)

# Compute cross-contrast concordance
if (!is.null(attrib_c1)) {
  # C1 vs C2 comparison
  c1_c2_genes <- intersect(attrib_c1$gene, attrib_c2$gene)
  cat(sprintf("\n  C1 vs C2 overlap: %d genes\n", length(c1_c2_genes)))

  c1_cats <- attrib_c1[gene %in% c1_c2_genes, setNames(category, gene)]
  c2_cats <- attrib_c2[gene %in% c1_c2_genes, setNames(category, gene)]

  # Genes intrinsic in C1 but composition-driven in C2 (or vice versa)
  c1_intrinsic_c2_comp <- sum(c1_cats[c1_c2_genes] == "Hepatocyte_intrinsic" &
                              c2_cats[c1_c2_genes] == "Composition_driven", na.rm = TRUE)
  c1_comp_c2_intrinsic <- sum(c1_cats[c1_c2_genes] == "Composition_driven" &
                              c2_cats[c1_c2_genes] == "Hepatocyte_intrinsic", na.rm = TRUE)
  both_intrinsic <- sum(c1_cats[c1_c2_genes] == "Hepatocyte_intrinsic" &
                        c2_cats[c1_c2_genes] == "Hepatocyte_intrinsic", na.rm = TRUE)

  cat(sprintf("  Both intrinsic: %d\n", both_intrinsic))
  cat(sprintf("  C1 intrinsic -> C2 composition-driven: %d\n", c1_intrinsic_c2_comp))
  cat(sprintf("  C1 composition-driven -> C2 intrinsic: %d\n", c1_comp_c2_intrinsic))

  # Key question: does NAFL->NASH show MORE immune-driven signal than Disease-vs-Control?
  c1_degs <- attrib_c1[sig_unadj == TRUE]
  c2_degs <- attrib_c2[sig_unadj == TRUE]
  c1_comp_pct <- 100 * sum(c1_degs$category == "Composition_driven") / nrow(c1_degs)
  c2_comp_pct <- 100 * sum(c2_degs$category == "Composition_driven") / nrow(c2_degs)
  cat(sprintf("\n  Composition-driven fraction of DEGs:\n"))
  cat(sprintf("    C1 (Disease vs Control): %.1f%%\n", c1_comp_pct))
  cat(sprintf("    C2 (NASH vs NAFL):       %.1f%%\n", c2_comp_pct))
  if (c2_comp_pct > c1_comp_pct) {
    cat(sprintf("    --> NAFL->NASH shows %.1f%% MORE immune-driven signal than Disease-vs-Control\n",
      c2_comp_pct - c1_comp_pct))
  } else {
    cat(sprintf("    --> Disease-vs-Control shows %.1f%% MORE immune-driven signal than NAFL->NASH\n",
      c1_comp_pct - c2_comp_pct))
  }
}

# C2 vs C3 comparison
c2_c3_genes <- intersect(attrib_c2$gene, attrib_c3$gene)
cat(sprintf("\n  C2 vs C3 overlap: %d genes\n", length(c2_c3_genes)))

c2_cats <- attrib_c2[gene %in% c2_c3_genes, setNames(category, gene)]
c3_cats <- attrib_c3[gene %in% c2_c3_genes, setNames(category, gene)]

c2_intrinsic_c3_comp <- sum(c2_cats[c2_c3_genes] == "Hepatocyte_intrinsic" &
                            c3_cats[c2_c3_genes] == "Composition_driven", na.rm = TRUE)
c3_intrinsic_c2_comp <- sum(c3_cats[c2_c3_genes] == "Hepatocyte_intrinsic" &
                            c2_cats[c2_c3_genes] == "Composition_driven", na.rm = TRUE)

c2_degs_sub <- attrib_c2[sig_unadj == TRUE]
c3_degs_sub <- attrib_c3[sig_unadj == TRUE]
c2_comp_pct_sub <- 100 * sum(c2_degs_sub$category == "Composition_driven") / max(nrow(c2_degs_sub), 1)
c3_comp_pct_sub <- 100 * sum(c3_degs_sub$category == "Composition_driven") / max(nrow(c3_degs_sub), 1)

cat(sprintf("  Composition-driven fraction of DEGs:\n"))
cat(sprintf("    C2 (NASH vs NAFL):         %.1f%%\n", c2_comp_pct_sub))
cat(sprintf("    C3 (Advanced vs Early fib): %.1f%%\n", c3_comp_pct_sub))

# Save comparison
fwrite(comparison_summary, file.path(ODIR, "progression_deconv_comparison.csv"))
cat("\n  Saved: progression_deconv_comparison.csv\n")

# ============================================================
#  Final Summary
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  DECONVOLUTION ATTRIBUTION — PROGRESSION CONTRASTS COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\nOutputs saved to:", ODIR, "\n")
cat("  1. progression_deconv_attribution_c2.csv\n")
cat("  2. progression_deconv_attribution_c3.csv\n")
cat("  3. progression_deconv_adjusted_c2.csv\n")
cat("  4. progression_deconv_adjusted_c3.csv\n")
if (!is.null(inter_c2_result)) {
  cat("  5. progression_deconv_interaction_c2.csv\n")
} else {
  cat("  5. progression_deconv_interaction_c2.csv — SKIPPED (rank-deficient)\n")
}
cat("  6. progression_deconv_comparison.csv\n")

cat("\nC2 (NASH vs NAFL) attribution:\n")
print(attrib_c2[, .N, by = category][order(-N)])

cat("\nC3 (Advanced vs Early Fibrosis) attribution:\n")
print(attrib_c3[, .N, by = category][order(-N)])

cat("\nDone:", format(Sys.time()), "\n")
