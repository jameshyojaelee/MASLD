#!/usr/bin/env Rscript
# 03_integrate_counts.R
# ---------------------------------------------------------------------------
# Merge counts from all datasets, apply filterByExpr, RLE normalization.
# Output: results/integration/merged_dge.rds
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

# --- Load previously merged counts and metadata ---
merged <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))

# Ensure empty strings are NA
for (col in names(meta)) if (is.character(meta[[col]])) meta[get(col) == "", (col) := NA]

# Load QC report and filter to passing samples
qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
# Use pass_technical (PCA + libsize) instead of pass_all so that samples
# failing only the sex check are included in non-sex analyses (dream, meta,
# consensus, deconv). Sex-stratified analysis (script 26) applies its own
# pass_sex filter.
pass_ids <- qc[pass_technical == TRUE, sample_id]
merged <- merged[, pass_ids]
meta   <- meta[sample_id %in% pass_ids]

# CRITICAL: align meta row order to match merged column order (pass_ids order).
# Subsetting meta by %in% does NOT reorder rows, so positional assignment
# (dge$samples$x <- meta$x) would assign wrong labels without this step.
meta <- meta[match(pass_ids, meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(merged)))

cat("QC-passing samples:", ncol(merged), "\n")

# --- Create DGEList ---
dge <- DGEList(counts = merged)

# Add sample-level info
dge$samples$dataset  <- meta$dataset
dge$samples$condition <- factor(meta$condition)
dge$samples$group_binary <- factor(meta$group_binary, levels = c("Control", "Disease"))
# P1 fix 2026-05-28: store the RESOLVED sex (Script 01 `sex_final` = reported sex
# when available, k-means inferred otherwise), not the raw `meta$sex` which is NA
# for cohorts lacking annotated sex (GSE135251, GSE213621, GSE240729). Without
# this, any consumer reading dge$samples$sex from merged_dge.rds gets NAs.
# Prefer sex_final; fall back to inferred_sex, then raw sex.
dge$samples$sex <- if ("sex_final" %in% names(meta)) {
  meta$sex_final
} else if ("inferred_sex" %in% names(meta)) {
  meta$inferred_sex
} else {
  meta$sex
}

# --- Filter genes ---
# NOTE: filterByExpr includes all QC-passing samples (including disease-only
# cohorts) to maximize gene retention. The kallisto canonical pipeline uses
# tximport which handles gene universes differently.
design_filter <- model.matrix(~ 0 + group_binary, data = dge$samples)
keep <- filterByExpr(dge, design = design_filter)
cat("Genes passing filterByExpr:", sum(keep), "/", nrow(dge), "\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]

# --- RLE normalization ---
# RLE (DESeq2-style) is more robust than TMM when integrating studies with
# asymmetric group sizes. dream's (1|dataset) random effect handles
# inter-study differences, and RLE stabilizes normalization factors across
# unbalanced groups.
# (PRJNA512027 — with its L0/S0 library-prep confound — was permanently removed
# from the pipeline 2026-05-15.)
dge <- calcNormFactors(dge, method = "RLE")

# --- Summary ---
cat("\n===== MERGED DGE SUMMARY =====\n")
cat("Samples:", ncol(dge), "\n")
cat("Genes:", nrow(dge), "\n")
cat("Datasets:", paste(unique(dge$samples$dataset), collapse = ", "), "\n")
cat("Group distribution:\n")
print(table(dge$samples$group_binary, dge$samples$dataset))

# --- Save ---
saveRDS(dge, file.path(RDIR, "merged_dge.rds"))
cat("\nSaved:", file.path(RDIR, "merged_dge.rds"), "\n")
