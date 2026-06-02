#!/usr/bin/env bash
#SBATCH --job-name=compare_methods
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --partition=cpu
#SBATCH --output=RNA-seq/deconvolution/logs/compare_methods_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/compare_methods_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13

# Environment setup
MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv

echo "=== Running Deconvolution Methods Comparison ==="
"${MICROMAMBA}" run -p "${ENV_PREFIX}" \
  Rscript RNA-seq/deconvolution/scripts/16_compare_methods.R

echo ""
echo "=== Comparison Complete ==="
