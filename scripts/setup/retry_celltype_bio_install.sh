#!/bin/bash
# retry_celltype_bio_install.sh
# Targeted retry for failed Bioconductor + GitHub packages in celltype_bio env.
# Run after install_celltype_bio_env.sh. Uses longer timeouts + per-package
# error isolation.

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate celltype_bio

# Install cairo headers (needed by some plotting deps)
micromamba install -y -c conda-forge cairo pango harfbuzz fribidi libpng freetype || true

Rscript -e '
options(Ncpus = 4, timeout = 7200,
        repos = c(CRAN = "https://cloud.r-project.org"))

# Bioconductor
bioc_pkgs <- c("DESeq2","TOAST","speckle","MOFA2","SingleCellExperiment",
               "SummarizedExperiment","biomaRt","ComplexHeatmap",
               "OmnipathR","ComplexHeatmap")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("Installing Bioc: %s\n", p))
    res <- try(BiocManager::install(p, update = FALSE, ask = FALSE,
                                    force = FALSE, Ncpus = 4), silent = FALSE)
    if (inherits(res, "try-error")) cat("FAILED:", p, "\n")
  } else {
    cat(sprintf("%s: OK\n", p))
  }
}

# GitHub
gh_pkgs <- list(
  CARseq         = "Sun-lab/CARseq",
  MIND           = "randel/MIND",
  CellChat       = "jinworks/CellChat",
  nichenetr      = "saeyslab/nichenetr",
  multinichenetr = "saeyslab/multinichenetr",
  liana          = "saezlab/liana"
)
for (nm in names(gh_pkgs)) {
  if (!requireNamespace(nm, quietly = TRUE)) {
    cat(sprintf("Installing GitHub %s from %s...\n", nm, gh_pkgs[[nm]]))
    res <- try(devtools::install_github(gh_pkgs[[nm]], upgrade = "never",
                                        quiet = FALSE, force = FALSE), silent = FALSE)
    if (inherits(res, "try-error")) cat("FAILED:", nm, "\n")
  } else {
    cat(sprintf("%s: OK\n", nm))
  }
}

# Final status
cat("\n=== FINAL STATUS ===\n")
all_pkgs <- c("DESeq2","TOAST","speckle","MOFA2","OmnipathR",
              "CARseq","MIND","CellChat","nichenetr","multinichenetr","liana",
              "fgsea","edgeR","limma","decoupleR")
for (p in all_pkgs) cat(sprintf("%-18s %s\n", p,
                                ifelse(requireNamespace(p, quietly = TRUE), "OK", "MISSING")))
' 2>&1 | tail -100
