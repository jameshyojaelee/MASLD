#!/bin/bash
#SBATCH --job-name=fig_c13
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_c13_%j.log

module load R/4.4.0 2>/dev/null || true
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures
Rscript fig_c13_fibrosis_confounding.R
