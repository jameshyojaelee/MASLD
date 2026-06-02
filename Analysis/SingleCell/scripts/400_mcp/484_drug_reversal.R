#!/usr/bin/env Rscript
# 484_drug_reversal.R
# Per-program drug reversal analysis using existing LINCS/CGP results.

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
drug_dir <- file.path(root, "RNA-seq/results/drug_repurposing")

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]; k <- as.integer(args[2])
top_f <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_annot",
                   name, sprintf("program_topgenes.k%d.tsv", k))
out_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/validation/drug")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(top_f))

top <- fread(top_f)

# Candidate LINCS connectivity file
candidates <- c(
  file.path(drug_dir, "lincs_reversal_scores.tsv"),
  file.path(drug_dir, "cgp_reversal.tsv"),
  file.path(drug_dir, "lincs_fgsea.tsv"),
  Sys.glob(file.path(drug_dir, "*reversal*.tsv*"))
)
candidates <- candidates[file.exists(candidates)]
if (!length(candidates)) {
  cat("[484] No existing LINCS/CGP reversal files found; run Strategy 11/12 first\n")
  quit(status = 0)
}

results <- list()
for (f in candidates) {
  cat(sprintf("[484] reading %s\n", basename(f)))
  d <- fread(f)
  # fgsea-style: expect columns gene_set, pathway, NES, pval
  # If not applicable (score table), treat differently
  # Save quick preview
  print(head(d, 3))
}

# Per-program reversal requires running fgsea over LINCS signatures with the program gene set.
# Full implementation is a downstream step (Phase 8e); here we lay out inputs.
cat("[484] SKELETON -- full reversal test to be implemented after cNMF outputs materialize\n")
