#!/bin/bash
# Install Cicero + Monocle3 in the rnaseq R environment
# Run this interactively (NOT via sbatch) — needs network access

set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "Installing Monocle3 + Cicero in rnaseq env..."
echo "This may take 10-15 minutes."

Rscript -e '
# Install monocle3 dependencies first
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Monocle3 from GitHub (not yet on Bioconductor)
if (!requireNamespace("monocle3", quietly = TRUE)) {
    message("Installing monocle3...")
    BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment",
                           "DelayedArray", "DelayedMatrixStats",
                           "limma", "S4Vectors", "BiocGenerics",
                           "lme4", "terra", "ggrastr"),
                         update = FALSE, ask = FALSE)
    devtools::install_github("cole-trapnell-lab/monocle3", upgrade = "never")
}

# Cicero from GitHub
if (!requireNamespace("cicero", quietly = TRUE)) {
    message("Installing cicero...")
    devtools::install_github("cole-trapnell-lab/cicero-release", upgrade = "never")
}

# Verify
message("Checking installations...")
library(monocle3)
message("  monocle3: ", packageVersion("monocle3"))
library(cicero)
message("  cicero: ", packageVersion("cicero"))
message("Done!")
'
