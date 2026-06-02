#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --job-name=fig1_panels
#SBATCH --output=logs/fig1_integration_panels_%j.out
#SBATCH --error=logs/fig1_integration_panels_%j.err

set -euo pipefail

export PATH="/gpfs/commons/home/jameslee/miniforge3/condabin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/figures/fig1_integration_value_panels.R
