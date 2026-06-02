#!/bin/bash
#SBATCH --job-name=fig_pred_supp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

Rscript "${MASLD_PROJECT_ROOT}/scripts/figures/fig_prediction_supp.R"
