#!/usr/bin/env Rscript
# 02_per_study_de.R
# ---------------------------------------------------------------------------
# Per-study differential expression analysis using limma-voom.
# DE parameters (formula, contrast, condition_levels, norm_method) are read
# from config/human_datasets.yaml — no dataset-specific code blocks needed.
#
# TO ADD A NEW DATASET: add its de: section to config/human_datasets.yaml.
#
# Input:  merged_counts_raw.rds, meta_matched.rds, sample_qc_report.csv
# Output: results/per_study/{dataset}_de_results.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(yaml)
})

# Resolve project root
PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
BASE <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts")
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/per_study")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

# --- Load dataset config ---
cfg_path <- file.path(PROJECT_ROOT, "config/human_datasets.yaml")
if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
datasets_cfg <- yaml.load_file(cfg_path)$datasets

# --- Load data ---
merged <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
for (col in names(meta)) if (is.character(meta[[col]])) meta[get(col) == "", (col) := NA]
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))

# Filter to QC-passing samples
pass_ids <- qc[pass_technical == TRUE, sample_id]
cat("QC-passing samples:", length(pass_ids), "/", nrow(meta), "\n")
merged <- merged[, pass_ids]
meta   <- meta[sample_id %in% pass_ids]

# ============================================================================
# Core DE function — shared across all datasets
# ============================================================================
run_per_study_de <- function(dataset_name, design_formula, contrast_name,
                             contrast_levels, meta_sub, counts_sub,
                             norm_method = "TMM",
                             use_quality_weights = FALSE) {
  cat("\n", strrep("=", 60), "\n  ", dataset_name, "\n", strrep("=", 60), "\n")
  cat("  Samples:", ncol(counts_sub), "  Norm:", norm_method,
      "  QualityWeights:", use_quality_weights, "\n")

  # Create DGEList and filter
  dge <- DGEList(counts = counts_sub)
  design <- model.matrix(design_formula, data = meta_sub)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = norm_method)
  cat("  Genes after filterByExpr:", sum(keep), "\n")

  # Voom
  if (use_quality_weights) {
    v <- voomWithQualityWeights(dge, design, plot = FALSE)
  } else {
    v <- voom(dge, design, plot = FALSE)
  }

  # Fit
  fit <- lmFit(v, design)

  # Contrast
  if (!is.null(contrast_levels)) {
    contr <- makeContrasts(contrasts = contrast_name, levels = design)
    fit2 <- contrasts.fit(fit, contr)
  } else {
    fit2 <- fit
  }
  fit2 <- eBayes(fit2)

  # Extract results — add SE + df.total for downstream metafor/mashr arms
  # (mega_validation bundle, 2026-05-19).
  coef_idx <- if (is.null(contrast_levels)) contrast_name else 1
  res <- topTable(fit2, coef = coef_idx, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res$dataset <- dataset_name
  stdev_unsc <- if (is.matrix(fit2$stdev.unscaled)) fit2$stdev.unscaled[, coef_idx] else fit2$stdev.unscaled
  # Moderated SE: uses eBayes posterior variance (s2.post). Kept for reference but
  # must NOT be fed to metafor — moderated SE biases tau^2 downward in REML.
  res$SE_moderated   <- stdev_unsc * sqrt(fit2$s2.post)
  # Unmoderated SE: uses per-gene residual SD before eBayes shrinkage (sigma).
  # This is the correct SE for random-effects meta-analysis (Script 06).
  res$SE_unmoderated <- stdev_unsc * fit2$sigma
  # Keep legacy 'SE' column pointing to MODERATED for backward compat with any
  # downstream code that reads it, but new consumers should use SE_unmoderated.
  res$SE       <- res$SE_moderated
  res$df.total <- fit2$df.total
  res$df.residual <- fit2$df.residual
  res_dt <- as.data.table(res)

  # Summary
  sig <- res_dt[adj.P.Val < 0.1]
  cat("  DEGs (padj < 0.1):", nrow(sig),
      " (Up:", sig[logFC > 0, .N], " Down:", sig[logFC < 0, .N], ")\n")

  # Save
  outfile <- file.path(RDIR, paste0(dataset_name, "_de_results.csv"))
  fwrite(res_dt, outfile)
  cat("  Saved:", outfile, "\n")

  return(res_dt)
}

# Null-coalescing operator for safe config field access
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ============================================================================
# Config-driven dispatch — no dataset-specific blocks needed here.
# Each dataset's formula, contrast, and norm parameters come from config.
# ============================================================================
all_results <- list()

for (ds in names(datasets_cfg)) {
  de_cfg <- datasets_cfg[[ds]]$de

  msk <- meta$dataset == ds
  if (sum(msk) == 0) {
    cat("\n  SKIPPING", ds, "— no samples in merged counts (featureCounts pending?)\n")
    next
  }

  meta_sub <- as.data.frame(meta[msk])
  counts_sub <- merged[, meta[msk, sample_id]]

  # Set condition factor with dataset-specific levels from config
  meta_sub$condition <- factor(meta_sub$condition, levels = de_cfg$condition_levels)

  # Coerce available covariates
  # age: coerce to numeric if the column exists and has at least some non-NA values
  if ("age" %in% names(meta_sub) && !all(is.na(meta_sub$age))) {
    meta_sub$age <- as.numeric(meta_sub$age)
  }
  # sex: coerce to factor if the column exists and has at least some non-NA values
  if ("sex" %in% names(meta_sub) && !all(is.na(meta_sub$sex))) {
    meta_sub$sex <- factor(meta_sub$sex)
  }
  # inferred_sex: populate for datasets where sex_annotated=false (GSE135251, GSE213621).
  # For datasets with sex_annotated=true, inferred_sex is available as a validation column
  # but is not used in the DE formula (annotated sex is preferred there).
  # Note which source is used and warn the user.
  # NOTE: sex_annotated is at the top-level dataset config, NOT inside de:
  # Use datasets_cfg[[ds]]$sex_annotated, not de_cfg$sex_annotated
  sex_annotated_flag <- isTRUE(datasets_cfg[[ds]]$sex_annotated %||% TRUE)
  if (!sex_annotated_flag) {
    if ("inferred_sex" %in% names(meta_sub) && !all(is.na(meta_sub$inferred_sex))) {
      meta_sub$inferred_sex <- factor(meta_sub$inferred_sex)
      cat("  NOTE:", ds, "— sex NOT annotated; using XIST/DDX3Y k-means inferred_sex as covariate\n")
      cat("         inferred_sex distribution:", paste(table(meta_sub$inferred_sex), collapse = "/"), "(M/F)\n")
    } else {
      cat("  WARNING:", ds, "— sex_annotated=false but inferred_sex missing; dropping sex from formula\n")
    }
  }

  # Drop rows where required covariates are NA
  # (safety: formula references only columns that exist in meta_sub)
  formula_vars <- all.vars(as.formula(de_cfg$formula))
  formula_vars <- formula_vars[formula_vars != "0"]  # strip intercept zero
  for (v in formula_vars) {
    if (v %in% names(meta_sub)) {
      n_before <- nrow(meta_sub)
      meta_sub <- meta_sub[!is.na(meta_sub[[v]]), ]
      if (nrow(meta_sub) < n_before) {
        cat("  NOTE:", ds, "— dropped", n_before - nrow(meta_sub),
            "samples with NA in covariate:", v, "\n")
        counts_sub <- counts_sub[, meta_sub$sample_id]
      }
    }
  }

  tryCatch({
    res <- run_per_study_de(
      dataset_name       = ds,
      design_formula     = as.formula(de_cfg$formula),
      contrast_name      = de_cfg$contrast,
      contrast_levels    = TRUE,
      meta_sub           = meta_sub,
      counts_sub         = counts_sub,
      norm_method        = de_cfg$norm_method,
      use_quality_weights = de_cfg$quality_weights
    )
    # Annotate with sex source so downstream scripts can track provenance
    # Use top-level sex_annotated flag (not de_cfg, which doesn't have this field)
    res[, sex_source := ifelse(isTRUE(datasets_cfg[[ds]]$sex_annotated %||% TRUE),
                               "annotated", "inferred_kmeans")]
    # Overwrite saved file with annotated version
    fwrite(res, file.path(RDIR, paste0(ds, "_de_results.csv")))
    all_results[[ds]] <- res
  }, error = function(e) {
    cat("  ERROR in", ds, ":", conditionMessage(e), "\n")
  })
}

# ============================================================================
# Summary across datasets
# ============================================================================
cat("\n", strrep("=", 60), "\n")
cat("  CROSS-DATASET DE SUMMARY\n")
cat(strrep("=", 60), "\n\n")

summary_dt <- rbindlist(lapply(all_results, function(r) {
  data.table(
    dataset   = r$dataset[1],
    total_genes = nrow(r),
    degs_01   = sum(r$adj.P.Val < 0.1),
    degs_up   = sum(r$adj.P.Val < 0.1 & r$logFC > 0),
    degs_down = sum(r$adj.P.Val < 0.1 & r$logFC < 0)
  )
}))
print(summary_dt)

fwrite(summary_dt, file.path(RDIR, "per_study_summary.csv"))
cat("\nPer-study DE complete.\n")
