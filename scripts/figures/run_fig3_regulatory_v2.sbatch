#!/bin/bash
#SBATCH --job-name=fig3_reg_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3_reg_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3_reg_v2_%j.err

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

mkdir -p scripts/figures/logs

echo "[fig3_v2] start $(date)"
Rscript scripts/figures/fig3_regulatory_architecture_v2.R
echo "[fig3_v2] done  $(date)"
