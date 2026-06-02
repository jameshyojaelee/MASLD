#!/bin/bash
#SBATCH --job-name=fig_prediction
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_prediction_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_prediction_%j.err
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

set -euo pipefail

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${MASLD_PROJECT_ROOT}"

echo "=== Starting fig_prediction_models.R ==="
echo "Date: $(date)"
echo "Node: $(hostname)"

Rscript scripts/figures/fig_prediction_models.R

echo "=== Done ==="
echo "Date: $(date)"
