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

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
BASE <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts")
SOURCE_INT <- file.path(BASE, "analysis/integration")
RUN_ROOT <- Sys.getenv("MASLD_RUN_ROOT", "")
REMEDIATION_MODE <- nzchar(RUN_ROOT)
if (REMEDIATION_MODE) {
  RUN_ROOT <- normalizePath(RUN_ROOT, mustWork = TRUE)
  if (!file.exists(file.path(RUN_ROOT, ".bg001_candidate_root"))) stop("Missing BG-001 candidate sentinel")
  INPUT_ROOT <- Sys.getenv("MASLD_INPUT_RUN_ROOT", RUN_ROOT)
  INPUT_ROOT <- normalizePath(INPUT_ROOT, mustWork = TRUE)
  if (!file.exists(file.path(INPUT_ROOT, ".bg001_candidate_root"))) stop("Input run root lacks BG-001 sentinel")
  INT <- RUN_ROOT
  INPUT_INT <- INPUT_ROOT
} else {
  INT <- SOURCE_INT
  INPUT_INT <- SOURCE_INT
}
RDIR <- file.path(INT, "results/integration")
INPUT_RDIR <- file.path(INPUT_INT, "results/integration")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

# --- Load previously merged counts and metadata ---
merged <- readRDS(file.path(INPUT_RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(INPUT_RDIR, "meta_matched.rds"))

# Ensure empty strings are NA
for (col in names(meta)) if (is.character(meta[[col]])) meta[get(col) == "", (col) := NA]

# Load QC report and filter to passing samples
qc <- fread(file.path(INPUT_INT, "qc/sample_qc_report.csv"))
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

filter_scope <- Sys.getenv("MASLD_FILTER_SCOPE", "legacy_all")
if (!filter_scope %in% c("legacy_all", "canonical_five")) {
  stop("MASLD_FILTER_SCOPE must be legacy_all or canonical_five")
}
counts_mode <- Sys.getenv("MASLD_COUNTS_MODE", if (REMEDIATION_MODE) "" else "current")
qc_mode <- Sys.getenv("MASLD_QC_MODE", if (REMEDIATION_MODE) "" else "current")
if (REMEDIATION_MODE && !counts_mode %in% c("read", "fragment")) {
  stop("MASLD_COUNTS_MODE must be read or fragment in remediation mode")
}
if (REMEDIATION_MODE && !qc_mode %in% c("recompute", "locked", "reuse_F_legacy")) {
  stop("MASLD_QC_MODE does not name a BG-001 remediation QC mode")
}
if (filter_scope == "canonical_five") {
  cfg_path <- Sys.getenv("MASLD_CONFIG_PATH", file.path(PROJECT_ROOT, "config/human_datasets.yaml"))
  cfg <- yaml::read_yaml(cfg_path)$datasets
  mega <- names(Filter(function(x) isTRUE(x$de$include_in_mega), cfg))
  keep_samples <- meta$dataset %in% mega
  merged <- merged[, keep_samples, drop = FALSE]
  meta <- meta[keep_samples]
  if (length(unique(meta$dataset)) != 5L) stop("Canonical-five scope did not resolve exactly five cohorts")
  stopifnot(all(meta$sample_id == colnames(merged)))
}

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
gene_mode <- Sys.getenv("MASLD_GENE_MODE", "native")
if (!gene_mode %in% c("native", "locked")) stop("MASLD_GENE_MODE must be native or locked")
if (REMEDIATION_MODE) {
  raw_logcpm_mean <- rowMeans(cpm(dge, log = TRUE, prior.count = 1))
  cohort_support <- integer(nrow(dge))
  for (dataset_id in unique(dge$samples$dataset)) {
    cohort_support <- cohort_support + as.integer(rowSums(dge$counts[, dge$samples$dataset == dataset_id, drop = FALSE]) > 0)
  }
  filter_statistics <- data.table(
    gene = rownames(dge),
    total_count = rowSums(dge$counts),
    n_samples_nonzero = rowSums(dge$counts > 0),
    n_cohorts_nonzero = cohort_support,
    mean_logCPM = raw_logcpm_mean
  )
}
if (gene_mode == "locked") {
  locked_dge_path <- Sys.getenv("MASLD_LOCKED_DGE", file.path(SOURCE_INT, "results/integration/merged_dge.rds"))
  locked_genes <- rownames(readRDS(locked_dge_path))
  missing_genes <- setdiff(locked_genes, rownames(dge))
  if (length(missing_genes)) stop("Candidate matrix lacks ", length(missing_genes), " locked genes")
  dge <- dge[locked_genes, , keep.lib.sizes = FALSE]
  cat("Genes locked to canonical universe:", nrow(dge), "\n")
} else {
  keep <- filterByExpr(dge, design = design_filter)
  cat("Genes passing filterByExpr:", sum(keep), "/", nrow(dge), "\n")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
}
if (REMEDIATION_MODE) {
  filter_statistics[, kept := gene %in% rownames(dge)]
  fwrite(
    filter_statistics,
    file.path(INT, "provenance/gene_filter_statistics.tsv.gz"),
    sep = "\t",
    compress = "gzip"
  )
}

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

if (REMEDIATION_MODE) {
  provenance_dir <- file.path(INT, "provenance")
  dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
  cohort_counts <- as.data.table(as.data.frame.matrix(table(dge$samples$dataset, dge$samples$group_binary)), keep.rownames = "dataset")
  fwrite(cohort_counts, file.path(provenance_dir, "cohort_group_counts.tsv"), sep = "\t")
  normalization <- as.data.table(dge$samples, keep.rownames = "sample_id")
  normalization[, effective_library_size := lib.size * norm.factors]
  fwrite(normalization, file.path(provenance_dir, "normalization_factors.tsv"), sep = "\t")
  manifest <- data.table(
    counts_mode = counts_mode,
    qc_mode = qc_mode,
    filter_scope = filter_scope,
    gene_mode = gene_mode,
    n_samples = ncol(dge),
    n_genes = nrow(dge),
    n_cohorts = length(unique(dge$samples$dataset)),
    norm_factor_min = min(dge$samples$norm.factors),
    norm_factor_median = median(dge$samples$norm.factors),
    norm_factor_max = max(dge$samples$norm.factors)
  )
  fwrite(manifest, file.path(provenance_dir, "preprocessing_manifest.tsv"), sep = "\t")
  writeLines(capture.output(sessionInfo()), file.path(provenance_dir, "integration_session_info.txt"))
}
