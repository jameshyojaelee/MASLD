#!/bin/bash
# install_celltype_bio_env.sh
#
# Creates conda env `celltype_bio` with Tier 1 dependencies for the 27-analysis
# cell-type-resolved MASLD biology pipeline (spec: docs/superpowers/specs/
# 2026-04-17-cell-type-resolved-masld-biology-design.md).
#
# Run from an interactive CPU node — NOT the login node.
#   srun --partition=cpu --mem=16G --cpus-per-task=4 --time=4:00:00 --pty bash
#   bash scripts/setup/install_celltype_bio_env.sh
#
# Packages: CARseq, bMIND, TOAST, speckle(propeller), CellChat v2, nichenetr,
# multinichenetr, MOFA2, liana(R), OmnipathR, decoupleR, fgsea (reuse).
#
# CellPhoneDB v5, Tensor-cell2cell, scCODA stay in `rapids_singlecell` (GPU env).
# COMMOT/stLearn/MISTy stay in `spatial`.

set -euo pipefail

ENV_NAME="celltype_bio"
R_VERSION="4.4"

echo "[1/4] Creating conda env $ENV_NAME (R $R_VERSION)..."
# Strategy: minimal base env (R 4.4 only); install Bioc/GitHub packages via R
# afterwards. Avoids conda dep-resolution conflicts with
# bioconductor-omnipathr 3.18.4 which requires R >= 4.5.
micromamba create -y -n "$ENV_NAME" -c conda-forge \
  r-base=${R_VERSION} \
  r-devtools r-biocmanager r-remotes r-data.table r-ggplot2 r-patchwork \
  r-matrix r-dplyr r-tidyr r-igraph r-nmf r-rcolorbrewer r-circlize \
  r-seurat r-seuratobject \
  r-optparse r-yaml r-jsonlite r-curl r-openssl \
  cmake pkg-config || echo "WARN: partial install, continuing..."

echo "[2/4] Activating env and installing Bioconductor + GitHub packages..."
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate "$ENV_NAME"

Rscript -e '
options(Ncpus = 4, timeout = 1200, repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor packages (install first; GitHub packages may depend on them)
bioc_pkgs <- c("limma","edgeR","DESeq2","fgsea","OmnipathR","decoupleR",
               "TOAST","speckle","MOFA2","SingleCellExperiment",
               "SummarizedExperiment","biomaRt","ComplexHeatmap")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("BiocManager::install %s\n", p))
    try(BiocManager::install(p, update = FALSE, ask = FALSE))
  }
}

# GitHub packages
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
    cat(sprintf("Installing %s from %s...\n", nm, gh_pkgs[[nm]]))
    try(devtools::install_github(gh_pkgs[[nm]], upgrade = "never", quiet = TRUE))
  } else {
    cat(sprintf("%-18s already installed\n", nm))
  }
}
' 2>&1 | tail -80

echo "[3/4] Verifying Tier 1 package availability..."
Rscript -e '
pkgs <- c("CARseq","MIND","CellChat","nichenetr","multinichenetr","liana",
          "TOAST","speckle","MOFA2","fgsea","edgeR","limma","DESeq2",
          "OmnipathR","decoupleR")
status <- sapply(pkgs, requireNamespace, quietly = TRUE)
cat(sprintf("%-18s %s\n", names(status), ifelse(status, "OK", "MISSING")))
if (any(!status)) {
  cat("\nMISSING packages:", paste(names(status)[!status], collapse=", "), "\n")
  cat("Install manually via BiocManager::install() or devtools::install_github().\n")
}'

echo "[4/4] Done. Activate env: micromamba activate $ENV_NAME"
