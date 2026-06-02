#!/bin/bash
#SBATCH --job-name=fig_pred_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_prediction_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_prediction_v2_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
mkdir -p "${BASE}/scripts/figures/logs"

source ~/.bashrc
micromamba activate rnaseq

cd "${BASE}"
Rscript scripts/figures/fig_prediction_v2.R
