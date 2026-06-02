#!/bin/bash
#SBATCH --job-name=fig
#SBATCH --partition=test
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/%x_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Started: $(date)"
echo "Running: scripts/figures/${RSCRIPT}"
Rscript "scripts/figures/${RSCRIPT}"
echo "Finished: $(date)"
