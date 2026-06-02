#!/usr/bin/env Rscript
# dream_loo_cv_contrasts.R
# ---------------------------------------------------------------------------
# Generalized leave-one-cohort-out cross-validation for three contrasts:
#
#   CONTRAST = "disease_vs_control"  Disease vs Control mega-analysis on the 5
#                                    yaml include_in_mega cohorts (GSE126848,
#                                    GSE130970, GSE135251, GSE162694, GSE213621).
#                                    Mirrors the original dream_loo_cv.R.
#
#   CONTRAST = "mash_vs_masl"        NASH (+ Borderline) vs NAFL on 7 cohorts
#                                    (GSE126848, GSE130970, GSE135251, GSE162694,
#                                    GSE167523, GSE174478, GSE193066).
#                                    MASH_DEF = "borderline_grouped" (default) or "strict".
#
#   CONTRAST = "mash_vs_healthy"     NASH (+ Borderline) vs Control on 4 cohorts
#                                    (GSE126848, GSE130970, GSE135251, GSE162694).
#                                    MASH_DEF as above.
#
# Env vars:
#   CONTRAST  required (one of the three above)
#   HELD_OUT  required (cohort id to hold out; must belong to the contrast pool)
#   MASH_DEF  optional, default "borderline_grouped". Only used for MASH contrasts.
#
# GENE-UNIVERSE-PRESERVATION INVARIANT:
#   Following 05_dream_mega_analysis.R / dream_loo_cv.R, calcNormFactors() and
#   filterByExpr() run on the FULL contrast cohort set BEFORE holding out a cohort.
#   Doing so keeps the gene universe and norm factors identical between the full
#   mega fit and every LOO fold so results are directly comparable.
#   voomWithDreamWeights() recomputes weights for the held-out subset.
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

contrast  <- Sys.getenv("CONTRAST", "")
held_out  <- Sys.getenv("HELD_OUT", "")
mash_def  <- Sys.getenv("MASH_DEF", "borderline_grouped")

if (nchar(contrast) == 0) {
  stop("CONTRAST env var must be set (disease_vs_control | mash_vs_masl | mash_vs_healthy | masl_vs_healthy)")
}
if (!contrast %in% c("disease_vs_control", "mash_vs_masl", "mash_vs_healthy", "masl_vs_healthy")) {
  stop("CONTRAST must be one of: disease_vs_control, mash_vs_masl, mash_vs_healthy, masl_vs_healthy")
}
if (nchar(held_out) == 0) stop("HELD_OUT env var must be set (e.g., GSE213621)")
if (!mash_def %in% c("borderline_grouped", "strict")) {
  stop("MASH_DEF must be one of: borderline_grouped, strict")
}

cat("=== Dream LOO-CV ===\n")
cat("CONTRAST:", contrast, "\n")
if (contrast %in% c("mash_vs_masl", "mash_vs_healthy")) cat("MASH_DEF:", mash_def, "\n")
cat("HELD_OUT:", held_out, "\n")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(yaml)
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

PROJECT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE    <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts")
INT     <- file.path(BASE, "analysis/integration")
RDIR    <- file.path(INT, "results/integration")

# --- Cohort pool per contrast ---
if (contrast == "disease_vs_control") {
  yaml_path <- file.path(PROJECT, "config/human_datasets.yaml")
  ycfg <- yaml::read_yaml(yaml_path)$datasets
  pool_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
} else if (contrast == "mash_vs_masl") {
  # PRJNA512027 permanently removed from pipeline 2026-05-15 (L0/S0 library
  # batch perfectly confounded with disease severity).
  pool_cohorts <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694",
                    "GSE167523", "GSE174478", "GSE193066")
} else {  # mash_vs_healthy OR masl_vs_healthy — 4 cohorts with Control + disease arm
  pool_cohorts <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694")
}

if (!held_out %in% pool_cohorts) {
  stop(held_out, " is not in the ", contrast, " cohort pool. Choose from: ",
       paste(pool_cohorts, collapse = ", "))
}

# --- Output directory ---
out_subdir <- contrast
if (contrast %in% c("mash_vs_masl", "mash_vs_healthy") && mash_def == "strict") {
  out_subdir <- paste0(contrast, "_strict")
}
loo_dir <- file.path(RDIR, "loo_cv", out_subdir)
dir.create(loo_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
#  Branch A: Disease vs Control (mirror dream_loo_cv.R, load merged_dge.rds)
# =============================================================================
if (contrast == "disease_vs_control") {

  cat("Loading merged DGE...\n")
  dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

  # NOTE: filterByExpr / calcNormFactors NOT re-run after subsetting — gene
  # universe and norm factors stay identical to the full mega fit.
  keep_samples <- dge$samples$dataset %in% setdiff(pool_cohorts, held_out)
  dge_loo <- dge[, keep_samples]

  meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
  stopifnot("All DGE samples must be in meta_matched" =
              all(colnames(dge_loo) %in% meta_new$sample_id))
  matched_sex <- meta_new$inferred_sex[match(colnames(dge_loo), meta_new$sample_id)]
  n_na_sex <- sum(is.na(matched_sex))
  if (n_na_sex > 0) {
    warning(n_na_sex, " samples have NA inferred_sex; these will be dropped by dream().")
  }

  info <- data.frame(
    group_binary = factor(dge_loo$samples$group_binary, levels = c("Control", "Disease")),
    dataset      = droplevels(factor(dge_loo$samples$dataset)),
    inferred_sex = factor(matched_sex),
    stringsAsFactors = FALSE
  )
  rownames(info) <- colnames(dge_loo)

  cat("Pool cohorts (k =", length(pool_cohorts), "):",
      paste(pool_cohorts, collapse = ", "), "\n")
  cat("Datasets remaining:", paste(sort(unique(as.character(info$dataset))), collapse = ", "), "\n")
  cat("Samples remaining:", ncol(dge_loo), "\n")
  cat("Genes (same as full model):", nrow(dge_loo), "\n")
  cat("\nGroup distribution:\n")
  print(table(info$group_binary, info$dataset))
  cat("Sex distribution:\n")
  print(table(info$inferred_sex, useNA = "always"))

  if (length(unique(info$dataset)) < 2) {
    stop("Need >= 2 datasets for random effect")
  }

  form <- ~ group_binary + inferred_sex + (1 | dataset)
  coef_name <- "group_binaryDisease"
  dge_final <- dge_loo
  info_final <- info

# =============================================================================
#  Branch B: MASH contrasts (build fresh DGE from merged_counts_raw.rds)
# =============================================================================
} else {

  cat("Loading raw counts + metadata...\n")
  counts <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
  meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
  qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
  meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

  # Restrict to the contrast's cohort pool
  meta <- meta[dataset %in% pool_cohorts]

  # Contrast-specific filter + factor construction
  if (contrast == "mash_vs_masl") {
    if (mash_def == "borderline_grouped") {
      meta_c <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline")]
      meta_c[, contrast_grp := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
    } else {
      meta_c <- meta[diagnosis_harmonized %in% c("NAFL", "NASH")]
      meta_c[, contrast_grp := as.character(diagnosis_harmonized)]
    }
    meta_c[, contrast_grp := factor(contrast_grp, levels = c("NAFL", "NASH"))]
    coef_name <- "contrast_grpNASH"
  } else if (contrast == "mash_vs_healthy") {
    if (mash_def == "borderline_grouped") {
      meta_c <- meta[diagnosis_harmonized %in% c("Control", "NASH", "Borderline")]
      meta_c[, contrast_grp := fifelse(diagnosis_harmonized == "Control", "Control", "MASH")]
    } else {
      meta_c <- meta[diagnosis_harmonized %in% c("Control", "NASH")]
      meta_c[, contrast_grp := fifelse(diagnosis_harmonized == "Control", "Control", "MASH")]
    }
    meta_c[, contrast_grp := factor(contrast_grp, levels = c("Control", "MASH"))]
    coef_name <- "contrast_grpMASH"
  } else {  # masl_vs_healthy — no Borderline ambiguity, single mode
    meta_c <- meta[diagnosis_harmonized %in% c("Control", "NAFL")]
    meta_c[, contrast_grp := factor(fifelse(diagnosis_harmonized == "Control", "Control", "MASL"),
                                    levels = c("Control", "MASL"))]
    coef_name <- "contrast_grpMASL"
  }

  cat("Pool cohorts (k =", length(pool_cohorts), "):",
      paste(pool_cohorts, collapse = ", "), "\n")
  cat("Contrast-eligible samples (full pool):", nrow(meta_c), "\n")
  cat("  Per dataset x group:\n")
  print(meta_c[, .N, by = .(dataset, contrast_grp)][order(dataset, contrast_grp)])

  # Build DGEList over the FULL contrast pool first
  meta_c <- meta_c[order(sample_id)]
  idx_all <- colnames(counts) %in% meta_c$sample_id
  dge_full <- DGEList(counts = counts[, idx_all])
  dge_full$samples <- cbind(dge_full$samples,
    meta_c[match(colnames(dge_full), meta_c$sample_id),
           .(contrast_grp, dataset, sex, inferred_sex, age)])

  # dataset_subbatch == dataset (PRJNA512027 L0/S0 subbatch logic permanently
  # removed 2026-05-15 with the cohort itself).
  dge_full$samples$dataset_subbatch <- factor(as.character(dge_full$samples$dataset))

  # sex_for_model: prefer annotated sex, fall back to inferred_sex
  dge_full$samples$sex_for_model <- dge_full$samples$sex
  na_sex <- is.na(dge_full$samples$sex_for_model) | dge_full$samples$sex_for_model == ""
  dge_full$samples$sex_for_model[na_sex] <- as.character(dge_full$samples$inferred_sex[na_sex])
  dge_full$samples$sex_for_model <- factor(dge_full$samples$sex_for_model)

  # Normalize + filter on the FULL contrast pool (gene-universe invariant)
  dge_full <- calcNormFactors(dge_full, method = "TMM")
  keep_full <- filterByExpr(dge_full, group = dge_full$samples$contrast_grp)
  dge_full <- dge_full[keep_full, , keep.lib.sizes = FALSE]
  cat("Genes after filterByExpr on full pool:", nrow(dge_full), "\n")
  cat("Samples in full pool:", ncol(dge_full), "\n")

  # NOW hold out the cohort — no re-normalize / re-filter
  keep_samples <- dge_full$samples$dataset != held_out
  dge_final <- dge_full[, keep_samples]
  dge_final$samples$dataset_subbatch <- droplevels(factor(dge_final$samples$dataset_subbatch))
  dge_final$samples$sex_for_model <- droplevels(factor(dge_final$samples$sex_for_model))
  dge_final$samples$contrast_grp <- droplevels(factor(
    dge_final$samples$contrast_grp,
    levels = levels(dge_full$samples$contrast_grp)))

  info_final <- data.frame(
    contrast_grp     = dge_final$samples$contrast_grp,
    dataset_subbatch = dge_final$samples$dataset_subbatch,
    sex_for_model    = dge_final$samples$sex_for_model,
    stringsAsFactors = FALSE
  )
  rownames(info_final) <- colnames(dge_final)

  cat("\nSamples after holding out", held_out, ":", ncol(dge_final), "\n")
  cat("Datasets remaining:",
      paste(sort(unique(as.character(dge_final$samples$dataset))), collapse = ", "), "\n")
  cat("\nGroup distribution (LOO):\n")
  print(table(info_final$contrast_grp, dge_final$samples$dataset))
  cat("Batch levels:",
      paste(levels(info_final$dataset_subbatch), collapse = ", "), "\n")
  cat("Sex levels:",
      paste(levels(info_final$sex_for_model), collapse = ", "), "\n")

  if (length(unique(info_final$dataset_subbatch)) < 2) {
    stop("Need >= 2 dataset_subbatch levels for random effect")
  }

  form <- ~ contrast_grp + sex_for_model + (1 | dataset_subbatch)
}

cat("\nFormula:", deparse(form), "\n")
cat("Coef:", coef_name, "\n")

# --- Parallelization ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

cat("Running voomWithDreamWeights...\n")
v <- suppressWarnings(voomWithDreamWeights(dge_final, form, info_final, BPPARAM = param))

cat("Running dream()...\n")
fit <- suppressWarnings(dream(v, form, info_final, BPPARAM = param))
# NOTE: do NOT call eBayes() after dream()

res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
res$gene <- rownames(res)
res_dt <- as.data.table(res)
setnames(res_dt, "adj.P.Val", "padj")

sig <- res_dt[padj < 0.1]
cat("\n===== LOO-CV RESULTS (", contrast, ", excluding ", held_out, ") =====\n", sep = "")
cat("Total genes tested:", nrow(res_dt), "\n")
cat("DEGs (padj < 0.1):", nrow(sig), "\n")
cat("  Up:", sum(sig$logFC > 0), "\n")
cat("  Down:", sum(sig$logFC < 0), "\n")

out_file <- file.path(loo_dir, paste0("dream_loo_", held_out, ".csv"))
fwrite(res_dt, out_file)
cat("\nSaved:", out_file, "\n")
cat("Done!\n")
