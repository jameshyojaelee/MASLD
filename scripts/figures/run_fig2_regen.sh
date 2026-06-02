#!/bin/bash
#SBATCH --job-name=fig2_regen
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig2_regen_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig2_regen_%j.log

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "[$(date)] starting fig2_progression_sex.R (Team F1 multi-step cascade pivot)"
Rscript scripts/figures/fig2_progression_sex.R
echo "[$(date)] done"
