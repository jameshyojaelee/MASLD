#!/usr/bin/env Rscript
# Convert 5-cohort canonical bulk RDS outputs to CSV for Python consumption.
#
# Source: RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/
# Files:
#   merged_counts_raw.rds   -> bulk_counts_raw.csv      (gene x sample, raw counts)
#   corrected_logcpm.rds    -> bulk_logcpm.csv           (gene x sample, log-CPM corrected)
#   merged_dge.rds          -> bulk_dge_metadata.csv     (sample-level metadata extracted from DGEList)
#
# Idempotent: skips files that already exist.
# Run via SLURM (NOT login node — files are 100MB+ in memory).

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SRC_DIR <- file.path(
  PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
)
DEST_DIR <- file.path(
  PROJECT_ROOT,
  "Analysis/Perturbation/data/pseudobulks"
)
dir.create(DEST_DIR, showWarnings = FALSE, recursive = TRUE)

convert <- function(src_name, dest_name, extractor) {
  src <- file.path(SRC_DIR, src_name)
  dest <- file.path(DEST_DIR, dest_name)
  if (file.exists(dest)) {
    cat(sprintf("[skip] %s already exists\n", dest))
    return(invisible())
  }
  if (!file.exists(src)) {
    cat(sprintf("[skip] source missing: %s\n", src))
    return(invisible())
  }
  cat(sprintf("[+] %s -> %s\n", src_name, dest_name))
  obj <- readRDS(src)
  df <- extractor(obj)
  fwrite(df, dest, row.names = TRUE)
  cat(sprintf("    wrote %d rows x %d cols\n", nrow(df), ncol(df)))
}

# Raw counts: assume matrix gene x sample
convert("merged_counts_raw.rds", "bulk_counts_raw.csv", function(obj) {
  if (is.matrix(obj)) {
    df <- as.data.frame(obj)
    df$gene <- rownames(obj)
    df <- df[, c("gene", setdiff(colnames(df), "gene"))]
    return(df)
  }
  # Some pipelines wrap counts in a DGEList
  if (inherits(obj, "DGEList")) {
    m <- obj$counts
    df <- as.data.frame(m)
    df$gene <- rownames(m)
    df <- df[, c("gene", setdiff(colnames(df), "gene"))]
    return(df)
  }
  stop("Unknown structure for merged_counts_raw.rds: ", class(obj)[1])
})

# Log-CPM: matrix gene x sample (already batch-corrected)
convert("corrected_logcpm.rds", "bulk_logcpm.csv", function(obj) {
  if (is.matrix(obj)) {
    df <- as.data.frame(obj)
    df$gene <- rownames(obj)
    df <- df[, c("gene", setdiff(colnames(df), "gene"))]
    return(df)
  }
  stop("Unknown structure for corrected_logcpm.rds: ", class(obj)[1])
})

# DGEList metadata: extract $samples
convert("merged_dge.rds", "bulk_dge_metadata.csv", function(obj) {
  if (inherits(obj, "DGEList")) {
    df <- as.data.frame(obj$samples)
    df$sample_id <- rownames(df)
    df <- df[, c("sample_id", setdiff(colnames(df), "sample_id"))]
    return(df)
  }
  stop("Unknown structure for merged_dge.rds: ", class(obj)[1])
})

cat("[done] convert_bulk_rds_to_csv complete\n")
