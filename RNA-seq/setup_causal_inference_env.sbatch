#!/bin/bash
#SBATCH --job-name=setup_causal_env
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/setup_causal_env_%j.out
#SBATCH --error=RNA-seq/logs/setup_causal_env_%j.err

# =============================================================================
# Phase 0C: Create Python environment for OTTERS / scDRS / FOCUS / CELLECT
# =============================================================================

set -uo pipefail  # No -e: pip warnings return nonzero, MAGMA download may fail
echo "=== Setting up causal_inference Python environment ==="
echo "Start: $(date)"

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"

# Create environment if it doesn't exist
if ! micromamba env list | grep -q "causal_inference"; then
  echo "Creating causal_inference environment..."
  micromamba create -n causal_inference python=3.10 -y -c conda-forge
else
  echo "causal_inference environment already exists"
fi

micromamba activate causal_inference

# Core scientific stack
pip install numpy scipy pandas scikit-learn matplotlib seaborn

# scDRS (single-cell disease relevance scoring)
echo "Installing scDRS..."
pip install scdrs

# MAGMA (download binary for gene-level GWAS analysis)
MAGMA_DIR="${BASE}/tools/magma"
if [[ ! -f "${MAGMA_DIR}/magma" ]]; then
  echo "Downloading MAGMA..."
  mkdir -p "${MAGMA_DIR}"
  cd "${MAGMA_DIR}"
  wget --retry-connrefused --tries=3 -q https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_static.zip || \
    wget --retry-connrefused --tries=3 -q https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip
  # Verify download
  if ! unzip -t magma_v1.10_static.zip >/dev/null 2>&1 && ! unzip -t magma_v1.10.zip >/dev/null 2>&1; then
    echo "MAGMA download corrupted. Skipping MAGMA — install manually later."
    cd "${BASE}"
    # Don't fail the whole script for MAGMA
  else
    unzip -o magma_v1.10_static.zip 2>/dev/null || unzip -o magma_v1.10.zip 2>/dev/null
  fi
  chmod +x magma
  cd "${BASE}"
  echo "MAGMA installed at ${MAGMA_DIR}/magma"
else
  echo "MAGMA already installed"
fi

# MAGMA reference data (1000G EUR, hg19)
MAGMA_REF="${MAGMA_DIR}/ref"
if [[ ! -f "${MAGMA_REF}/g1000_eur.bed" ]]; then
  echo "Downloading MAGMA 1000G EUR reference..."
  mkdir -p "${MAGMA_REF}"
  cd "${MAGMA_REF}"
  wget -q https://ctg.cncr.nl/software/MAGMA/ref/g1000_eur.zip
  unzip -o g1000_eur.zip
  cd "${BASE}"
  echo "MAGMA reference downloaded"
else
  echo "MAGMA reference already exists"
fi

# MAGMA gene annotation (hg19)
if [[ ! -f "${MAGMA_REF}/NCBI37.3.gene.loc" ]]; then
  echo "Downloading MAGMA gene locations..."
  cd "${MAGMA_REF}"
  wget -q https://ctg.cncr.nl/software/MAGMA/aux/NCBI37.3.gene.loc.gz
  gunzip NCBI37.3.gene.loc.gz
  cd "${BASE}"
fi

# CELLEX (for CELLECT expression specificity)
echo "Installing CELLEX..."
pip install cellex

# pybedtools (used by multiple tools)
pip install pybedtools

# anndata + scanpy (for scDRS h5ad handling)
pip install anndata scanpy

echo ""
echo "=== Verification ==="
python -c "
import scdrs; print(f'scDRS: {scdrs.__version__}')
import cellex; print(f'CELLEX: OK')
import numpy; print(f'numpy: {numpy.__version__}')
import scanpy; print(f'scanpy: {scanpy.__version__}')
print('All packages OK')
"

echo ""
echo "=== Environment setup complete ==="
echo "End: $(date)"
