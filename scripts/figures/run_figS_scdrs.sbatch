#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=4 --mem=96G --time=48:00:00
#SBATCH --job-name=figS_scdrs --output=scripts/figures/logs/figS_scdrs_%j.out --error=scripts/figures/logs/figS_scdrs_%j.err
set -eo pipefail; mkdir -p scripts/figures/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Step 1: Python preprocessing (parquet -> compact joined CSVs).
# Spatial env has pyarrow + pandas; rnaseq R env lacks arrow.
echo "[figS_scdrs] step 1: python prep"
micromamba activate spatial
python scripts/figures/figS_scdrs_prep.py
micromamba deactivate

# Step 2: R figure assembly (reads the prepped CSVs)
echo "[figS_scdrs] step 2: R figure assembly"
micromamba activate rnaseq
Rscript scripts/figures/figS_scdrs_bulk.R
Rscript scripts/figures/figS_scdrs_gwas.R

echo "[figS_scdrs] outputs:"
ls -la figures/supplementary/figS_scdrs/ || true
echo "[figS_scdrs] panels:"
ls figures/supplementary/figS_scdrs/panels/ || true
