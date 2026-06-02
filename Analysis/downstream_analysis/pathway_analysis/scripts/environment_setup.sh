#!/bin/bash
# Environment setup script for pathway analysis
# Creates a fresh conda environment with R and required packages

set -e

# Configuration
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
ENV_NAME="pathway_analysis"
ENV_PREFIX="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/.mamba/${ENV_NAME}"

echo "=== Pathway Analysis Environment Setup ==="
echo "Environment: ${ENV_PREFIX}"
echo ""

# Create environment with R and bioconda packages
echo "[1/4] Creating conda environment with R..."
$MICROMAMBA create -p "${ENV_PREFIX}" -c conda-forge -c bioconda -c r \
    r-base=4.4.3 \
    r-tidyverse \
    r-devtools \
    r-biocmanager \
    bioconductor-clusterprofiler \
    bioconductor-enrichplot \
    bioconductor-fgsea \
    bioconductor-reactomepa \
    bioconductor-dose \
    bioconductor-org.hs.eg.db \
    bioconductor-pathview \
    r-msigdbr \
    r-ggupset \
    r-ggridges \
    r-pheatmap \
    r-rcolorbrewer \
    r-scales \
    python=3.11 \
    pandas \
    matplotlib \
    seaborn \
    -y

echo ""
echo "[2/4] Verifying R package installations..."
$MICROMAMBA run -p "${ENV_PREFIX}" R --vanilla -e "
    pkgs <- c('clusterProfiler', 'enrichplot', 'fgsea', 'ReactomePA', 'DOSE', 'org.Hs.eg.db', 'pathview', 'msigdbr')
    for(pkg in pkgs) {
        if(require(pkg, character.only=TRUE)) {
            cat(paste0('[OK] ', pkg, '\n'))
        } else {
            cat(paste0('[MISSING] ', pkg, '\n'))
        }
    }
"

echo ""
echo "[3/4] Environment created successfully at: ${ENV_PREFIX}"
echo ""
echo "[4/4] To use this environment in scripts, add:"
echo "  MICROMAMBA=\"/gpfs/commons/home/jameslee/.local/bin/micromamba\""
echo "  ENV_PREFIX=\"${ENV_PREFIX}\""
echo "  \$MICROMAMBA run -p \"\${ENV_PREFIX}\" Rscript your_script.R"
echo ""
echo "=== Setup Complete ==="
