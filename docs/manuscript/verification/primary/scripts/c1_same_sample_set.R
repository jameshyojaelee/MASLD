#!/usr/bin/env Rscript
# Re-run C2 and C13 dream on the SAME sample set (615 samples).
# Purpose: separate the fibrosis-adjustment effect from the sample-set effect.
# If the 95.5% collapse mostly disappears when using the same sample set,
# then C1's headline number is driven by sample exclusion, NOT fibrosis adjustment.

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

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

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
ODIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts"
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

# --- Load data (same as scripts 13 and 130) ---
counts   <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta     <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc       <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta     <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

modeling_meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"))
meta <- merge(meta, modeling_meta[, .(sample_id, fibrosis_stage, nas_score,
  diagnosis_harmonized, nas_group, fib_ge3, nas_ge5)],
  by = "sample_id", all.x = TRUE, suffixes = c("", ".mm"))
for (col in c("fibrosis_stage", "nas_score", "diagnosis_harmonized")) {
  mm_col <- paste0(col, ".mm")
  if (mm_col %in% names(meta)) {
    na_mask <- is.na(meta[[col]]) | meta[[col]] == ""
    if (any(na_mask)) meta[[col]][na_mask] <- meta[[mm_col]][na_mask]
    meta[, (mm_col) := NULL]
  }
}
meta[, sex_covar := inferred_sex]
na_sex <- is.na(meta$sex_covar) | meta$sex_covar == ""
if (any(na_sex) && "sex" %in% names(meta)) {
  meta$sex_covar[na_sex] <- meta$sex[na_sex]
}
meta[, sex_covar := factor(sex_covar)]

# Subset to C13 samples
meta_c13 <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline") &
                 !is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta_c13[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
meta_c13[, nafl_nash := factor(nafl_nash, levels = c("NAFL", "NASH"))]
meta_c13[, fib_numeric := as.numeric(fibrosis_stage)]
ds_counts_c13 <- meta_c13[, .(n_nafl = sum(nafl_nash == "NAFL"),
                          n_nash = sum(nafl_nash == "NASH")), by = dataset]
valid_ds_c13 <- ds_counts_c13[n_nafl >= 2 & n_nash >= 2, dataset]
meta_c13 <- meta_c13[dataset %in% valid_ds_c13]
cat(sprintf("N = %d (NAFL %d, NASH %d)\n", nrow(meta_c13),
  sum(meta_c13$nafl_nash=="NAFL"), sum(meta_c13$nafl_nash=="NASH")))
cat(sprintf("Valid datasets: %s\n", paste(valid_ds_c13, collapse=", ")))

# Build DGE
keep_samples <- intersect(meta_c13$sample_id, colnames(counts))
meta_c13 <- meta_c13[sample_id %in% keep_samples]
dge <- DGEList(counts = counts[, keep_samples])
m <- meta_c13[match(colnames(dge), meta_c13$sample_id)]
dge$samples <- cbind(dge$samples, m[, .(dataset, sex_covar, nafl_nash, fib_numeric)])
dge$samples$dataset <- factor(dge$samples$dataset)
dge$samples$sex_covar <- factor(dge$samples$sex_covar)
dge$samples$nafl_nash <- factor(dge$samples$nafl_nash, levels = c("NAFL", "NASH"))
dge <- calcNormFactors(dge, method = "TMM")
keep_genes <- filterByExpr(dge, group = dge$samples$nafl_nash)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
cat(sprintf("Genes after filter: %d\n", nrow(dge)))

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
cat(sprintf("Using %d cores\n", ncpus))
BPPARAM <- if (ncpus > 1) MulticoreParam(ncpus, progressbar = FALSE) else SerialParam()

# Run dream WITHOUT fibrosis (C2-style on C13 sample set)
form_no_fib <- ~ nafl_nash + sex_covar + (1 | dataset)
cat("Running dream WITHOUT fibrosis on C13 sample set...\n")
t0 <- Sys.time()
vobj <- suppressWarnings(voomWithDreamWeights(dge, form_no_fib, dge$samples, BPPARAM = BPPARAM))
fit <- suppressWarnings(dream(vobj, form_no_fib, dge$samples, BPPARAM = BPPARAM))
tt <- topTable(fit, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
tt <- as.data.table(tt)
setnames(tt, "adj.P.Val", "padj")
cat(sprintf("  Done in %.1f min\n", as.numeric(Sys.time() - t0, units="mins")))

# Run dream WITH fibrosis (C13-style on C13 sample set)
form_fib <- ~ nafl_nash + fib_numeric + sex_covar + (1 | dataset)
cat("Running dream WITH fibrosis on C13 sample set...\n")
t0 <- Sys.time()
vobj2 <- suppressWarnings(voomWithDreamWeights(dge, form_fib, dge$samples, BPPARAM = BPPARAM))
fit2 <- suppressWarnings(dream(vobj2, form_fib, dge$samples, BPPARAM = BPPARAM))
tt2 <- topTable(fit2, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
tt2$gene <- rownames(tt2)
tt2 <- as.data.table(tt2)
setnames(tt2, "adj.P.Val", "padj")
cat(sprintf("  Done in %.1f min\n", as.numeric(Sys.time() - t0, units="mins")))

cat("\n=== COMPARISON: WITHOUT fib vs WITH fib, SAME sample set (615 samples) ===\n")
for (t in c(0.01, 0.05, 0.1)) {
  n_nofib <- sum(tt$padj < t, na.rm = TRUE)
  n_fib <- sum(tt2$padj < t, na.rm = TRUE)
  red <- 100 * (1 - n_fib/n_nofib)
  cat(sprintf("padj<%.2f: no_fib=%d, with_fib=%d, reduction=%.2f%%\n", t, n_nofib, n_fib, red))
}

# Save
fwrite(tt, file.path(ODIR, "c2_on_c13_samples.csv"))
fwrite(tt2, file.path(ODIR, "c13_reproduce.csv"))
cat(sprintf("\nSaved: %s/c2_on_c13_samples.csv\n", ODIR))
cat(sprintf("Saved: %s/c13_reproduce.csv\n", ODIR))

cat("\n=== Script complete ===\n")
