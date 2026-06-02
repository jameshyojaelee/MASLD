#!/bin/bash
#SBATCH --job-name=figS2_progression
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figS2_progression_%j.out
#SBATCH --error=logs/figS2_progression_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

Rscript scripts/figures/figS02_progression_contrasts.R
