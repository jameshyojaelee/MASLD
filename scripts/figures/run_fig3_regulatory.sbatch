#!/bin/bash
#SBATCH --job-name=fig3_deconf
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3_deconf_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3_deconf_%j.err

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
Rscript scripts/figures/fig3_deconfounding.R
