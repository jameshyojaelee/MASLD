#!/usr/bin/env Rscript
# V6-a Test 1: C1 via GeneMeta-style per-study limma + fibrosis adjustment
# =====================================================================
# Our claim: Fibrosis adjustment collapses 5,915 C2 DEGs to 268 C13 DEGs
# (95.5% reduction). Already verified at 90.16% on matched sample set
# (SLURM 15320822: 2,724 -> 268).
#
# Competitor (Piras & DiStefano 2024): random-effects meta-analysis
# (GeneMeta) + per-study limma with fibrosis as covariate within-study.
# Fits each dataset separately, then meta-analyzes with metafor (REML).
#
# Test: Does per-study limma NASH-vs-NAFL, adjusted for fibrosis
# within-study, meta-analyzed, also show >=90% DEG collapse from baseline?
#
# Uses same 615-sample C13 set (5 datasets: GSE130970, GSE135251,
# GSE162694, GSE174478, GSE193066) for direct comparability.
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(metafor)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUT_DIR <- file.path(BASE, "docs/manuscript/verification/method_comparison")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Load data ===\n")
counts <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

modeling_meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"))
meta <- merge(meta, modeling_meta[, .(sample_id, fibrosis_stage,
  diagnosis_harmonized)], by = "sample_id", all.x = TRUE,
  suffixes = c("", ".mm"))
for (col in c("fibrosis_stage", "diagnosis_harmonized")) {
  mm_col <- paste0(col, ".mm")
  if (mm_col %in% names(meta)) {
    na_mask <- is.na(meta[[col]]) | meta[[col]] == ""
    if (any(na_mask)) meta[[col]][na_mask] <- meta[[mm_col]][na_mask]
    meta[, (mm_col) := NULL]
  }
}
meta[, sex_covar := inferred_sex]
na_sex <- is.na(meta$sex_covar) | meta$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta)) meta$sex_covar[na_sex] <- meta$sex[na_sex]
meta[, sex_covar := factor(sex_covar)]

# Subset to C13-style sample set
meta_c13 <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline") &
                 !is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c13[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
meta_c13[, nafl_nash := factor(nafl_nash, levels = c("NAFL", "NASH"))]
meta_c13[, fib_numeric := as.numeric(fibrosis_stage)]
ds_counts <- meta_c13[, .(n_nafl = sum(nafl_nash == "NAFL"),
                         n_nash = sum(nafl_nash == "NASH")), by = dataset]
valid_ds <- ds_counts[n_nafl >= 2 & n_nash >= 2, dataset]
meta_c13 <- meta_c13[dataset %in% valid_ds]
cat(sprintf("N = %d (NAFL %d, NASH %d) across %d datasets\n",
  nrow(meta_c13), sum(meta_c13$nafl_nash=="NAFL"),
  sum(meta_c13$nafl_nash=="NASH"), length(valid_ds)))

# ============================================================================
# Run per-study limma (two flavors):
# (A) No fibrosis adjustment (baseline): y ~ nafl_nash + sex_covar
# (B) With fibrosis covariate:           y ~ nafl_nash + fib_numeric + sex_covar
# ============================================================================
run_per_study_limma <- function(ds, meta_ds, counts_ds, with_fib = FALSE) {
  # Build design matrix
  if (length(unique(meta_ds$sex_covar)) > 1) {
    if (with_fib) {
      design <- model.matrix(~ nafl_nash + fib_numeric + sex_covar, data = meta_ds)
    } else {
      design <- model.matrix(~ nafl_nash + sex_covar, data = meta_ds)
    }
  } else {
    if (with_fib) {
      design <- model.matrix(~ nafl_nash + fib_numeric, data = meta_ds)
    } else {
      design <- model.matrix(~ nafl_nash, data = meta_ds)
    }
  }

  # If fib covariate has <2 unique values, skip adjustment
  if (with_fib && length(unique(meta_ds$fib_numeric)) < 2) {
    cat(sprintf("  [%s] Only 1 fibrosis level; cannot adjust. Skipping.\n", ds))
    return(NULL)
  }

  dge <- DGEList(counts = counts_ds)
  dge <- calcNormFactors(dge, method = "TMM")
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 1000) {
    cat(sprintf("  [%s] Too few genes after filter (%d).\n", ds, nrow(dge)))
    return(NULL)
  }
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$dataset <- ds
  data.table(tt)
}

per_study_nofib <- list()
per_study_fib   <- list()
for (ds in valid_ds) {
  meta_ds <- meta_c13[dataset == ds]
  meta_ds <- as.data.frame(meta_ds)
  ids <- intersect(meta_ds$sample_id, colnames(counts))
  meta_ds <- meta_ds[match(ids, meta_ds$sample_id), ]
  counts_ds <- counts[, ids]
  cat(sprintf("\n--- %s (N=%d) ---\n", ds, nrow(meta_ds)))
  cat(sprintf("  NAFL=%d, NASH=%d, fib range: %s\n",
              sum(meta_ds$nafl_nash=="NAFL"),
              sum(meta_ds$nafl_nash=="NASH"),
              paste(range(meta_ds$fib_numeric), collapse="-")))
  tt0 <- run_per_study_limma(ds, meta_ds, counts_ds, with_fib = FALSE)
  tt1 <- run_per_study_limma(ds, meta_ds, counts_ds, with_fib = TRUE)
  if (!is.null(tt0)) per_study_nofib[[ds]] <- tt0
  if (!is.null(tt1)) per_study_fib[[ds]]   <- tt1
}

nofib_all <- rbindlist(per_study_nofib)
fib_all   <- rbindlist(per_study_fib)
cat(sprintf("\nPer-study (no fib): %d rows across %d datasets\n",
            nrow(nofib_all), length(per_study_nofib)))
cat(sprintf("Per-study (with fib): %d rows across %d datasets\n",
            nrow(fib_all), length(per_study_fib)))

# ============================================================================
# Random-effects meta-analysis (GeneMeta / metafor, REML)
# ============================================================================
meta_analyze <- function(per_study) {
  # Pivot to wide: one row per gene, columns per dataset
  logfc_wide <- dcast(per_study, gene ~ dataset, value.var = "logFC")
  se_wide    <- dcast(per_study, gene ~ dataset, value.var = "AveExpr") # placeholder
  # SE from t and logFC: SE = logFC / t
  per_study[, se := logFC / t]
  se_wide <- dcast(per_study, gene ~ dataset, value.var = "se")

  genes <- logfc_wide$gene
  logfc_mat <- as.matrix(logfc_wide[, -"gene"])
  se_mat    <- as.matrix(se_wide[, -"gene"])

  # Only meta-analyze genes present in >=3 datasets
  n_present <- rowSums(!is.na(logfc_mat) & !is.na(se_mat) & se_mat > 0)
  keep <- n_present >= 3
  cat(sprintf("  Genes in >=3 datasets: %d / %d\n", sum(keep), length(genes)))

  out <- data.table(gene = genes[keep])
  out[, c("meta_logFC", "meta_se", "meta_pval", "meta_Q", "meta_I2",
          "n_datasets") := NA_real_]

  logfc_keep <- logfc_mat[keep, , drop = FALSE]
  se_keep    <- se_mat[keep, , drop = FALSE]

  for (i in seq_len(nrow(out))) {
    yi <- logfc_keep[i, ]
    sei <- se_keep[i, ]
    valid <- !is.na(yi) & !is.na(sei) & sei > 0
    if (sum(valid) < 3) next
    res <- tryCatch(
      rma(yi = yi[valid], sei = sei[valid], method = "REML"),
      error = function(e) NULL
    )
    if (is.null(res)) next
    out[i, meta_logFC := res$b[1]]
    out[i, meta_se    := res$se]
    out[i, meta_pval  := res$pval]
    out[i, meta_Q     := res$QE]
    out[i, meta_I2    := res$I2]
    out[i, n_datasets := sum(valid)]
    if (i %% 2000 == 0) cat(sprintf("    gene %d / %d\n", i, nrow(out)))
  }

  out[, meta_padj := p.adjust(meta_pval, method = "BH")]
  out
}

cat("\n=== Meta-analysis (no fib adjustment) ===\n")
t0 <- Sys.time()
meta_nofib <- meta_analyze(nofib_all)
cat(sprintf("  %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
cat("\n=== Meta-analysis (with fib adjustment) ===\n")
t0 <- Sys.time()
meta_fib <- meta_analyze(fib_all)
cat(sprintf("  %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))

# Save
fwrite(meta_nofib, file.path(OUT_DIR, "test1_C1_genemeta_nofib.csv"))
fwrite(meta_fib,   file.path(OUT_DIR, "test1_C1_genemeta_fib.csv"))

# =====================================================================
# COMPARISON
# =====================================================================
cat("\n=== COMPARISON: GeneMeta per-study limma, NO fib vs WITH fib ===\n")
for (t in c(0.01, 0.05, 0.1)) {
  n_nofib <- sum(meta_nofib$meta_padj < t, na.rm = TRUE)
  n_fib <- sum(meta_fib$meta_padj < t, na.rm = TRUE)
  red <- 100 * (1 - n_fib/n_nofib)
  cat(sprintf("padj<%.2f: no_fib=%d, with_fib=%d, reduction=%.2f%%\n",
              t, n_nofib, n_fib, red))
}

# Our-method comparison
our_c2 <- 5915   # dream mega-analysis
our_c13 <- 268
our_same_nofib <- 2724 # from SLURM 15320822
our_same_fib <- 268
our_reduction_same <- 100*(1 - 268/2724)
cat(sprintf("\nOur method (dream, same 615 samples): %d -> %d = %.2f%% reduction\n",
            our_same_nofib, our_same_fib, our_reduction_same))

# Summary
n_nofib_05 <- sum(meta_nofib$meta_padj < 0.05, na.rm = TRUE)
n_fib_05 <- sum(meta_fib$meta_padj < 0.05, na.rm = TRUE)
comp_reduction <- 100 * (1 - n_fib_05/n_nofib_05)

summary_df <- data.table(
  method = c("dream (ours, same samples)", "GeneMeta per-study limma"),
  n_nofib_padj05 = c(2724, n_nofib_05),
  n_fib_padj05   = c(268, n_fib_05),
  pct_reduction  = c(our_reduction_same, comp_reduction),
  survives = c(our_reduction_same >= 80, comp_reduction >= 80)
)
print(summary_df)
fwrite(summary_df, file.path(OUT_DIR, "test1_C1_summary.csv"))

cat("\n=== Done ===\n")
