#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --job-name=auroc
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/cross_study_auroc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/cross_study_auroc_%j.err

# Cross-study AUROC prediction.
# Depends on all 8 dream LOO-CV jobs completing first.
# Usage: sbatch scripts/run_cross_study_auroc.sh

set -euo pipefail

export PATH="/gpfs/commons/home/jameslee/miniforge3/condabin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/cross_study_auroc.R
