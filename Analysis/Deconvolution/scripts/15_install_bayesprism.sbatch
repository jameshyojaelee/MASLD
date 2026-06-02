#!/usr/bin/env bash
#SBATCH --job-name=install_bayesprism
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=cpu
#SBATCH --output=RNA-seq/deconvolution/logs/install_bayesprism_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/install_bayesprism_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13

# Environment setup
export LD_LIBRARY_PATH="RNA-seq/.mamba/music_deconv/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export R_LIBS_USER=
export R_LIBS=
export R_LIBS_SITE="RNA-seq/.mamba/music_deconv/lib/R/library"

MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv

echo "=== Installing BayesPrism/InstaPrism packages ==="
"${MICROMAMBA}" run -p "${ENV_PREFIX}" \
  Rscript --vanilla RNA-seq/deconvolution/scripts/install_bayesprism.R

echo ""
echo "=== Installation Complete ==="
