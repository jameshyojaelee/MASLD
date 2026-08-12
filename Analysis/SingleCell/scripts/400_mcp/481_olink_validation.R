#!/usr/bin/env Rscript
# 481_olink_validation.R
# Validate cNMF programs in Olink plasma proteomics (1,461 proteins x 177 subjects).
# For each program's top secreted genes, test stage association in plasma.

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
olink_dir <- file.path(root, "data/olink")

# Locate Olink data (could be under multiple paths)
olink_candidates <- c(
  file.path(root, "data/olink/olink_data.rds"),
  file.path(root, "data/olink/olink_processed.tsv"),
  Sys.glob(file.path(root, "Analysis/Olink/results/olink_*.rds")),
  Sys.glob(file.path(root, "Analysis/Olink/results/olink_*.tsv.gz"))
)
olink_f <- olink_candidates[file.exists(olink_candidates)][1]
if (is.na(olink_f)) {
  cat("[481] No Olink data file found -- see docs/technical/DATASETS_AND_PIPELINES.md\n")
  quit(status = 0)
}
cat(sprintf("[481] Loading Olink: %s\n", olink_f))

# Accept either wide matrix (subjects x proteins) or long format
if (grepl("\\.rds$", olink_f)) {
  olink <- readRDS(olink_f)
} else {
  olink <- fread(olink_f)
}
print(str(olink))

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]; k <- as.integer(args[2])
top_f <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_annot",
                   name, sprintf("program_topgenes.k%d.tsv", k))
if (!file.exists(top_f)) {
  cat(sprintf("[481] cNMF annotation not found: %s\n", top_f))
  quit(status = 0)
}

top <- fread(top_f)
out_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/validation/olink")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Implementation placeholder: requires knowing the exact Olink file schema.
# This script is expected to be refined once the Olink file is located.
cat("[481] SKELETON script -- customize protein mapping once Olink schema is confirmed\n")
