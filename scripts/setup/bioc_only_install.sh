#!/bin/bash
# Focused Bioconductor-only install (no GitHub to avoid API rate limit).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate celltype_bio

Rscript -e '
options(Ncpus = 4, timeout = 7200,
        repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# One at a time with per-package error isolation
pkgs <- c("DESeq2","TOAST","speckle","MOFA2","biomaRt","ComplexHeatmap",
          "OmnipathR")
for (p in pkgs) {
  cat(sprintf("\n=== Installing %s ===\n", p))
  res <- try(BiocManager::install(p, update = FALSE, ask = FALSE, force = FALSE),
             silent = FALSE)
  status <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%s: %s\n", p, ifelse(status, "OK", "FAILED")))
}

cat("\n=== FINAL ===\n")
all_pkgs <- c("DESeq2","TOAST","speckle","MOFA2","OmnipathR","ComplexHeatmap",
              "fgsea","edgeR","limma","decoupleR")
for (p in all_pkgs) cat(sprintf("%-18s %s\n", p,
                                ifelse(requireNamespace(p, quietly = TRUE), "OK", "MISSING")))
' 2>&1
