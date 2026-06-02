#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --job-name=figS_B2
#SBATCH --output=scripts/figures/logs/figS_B2_liana_bulk_%j.out
#SBATCH --error=scripts/figures/logs/figS_B2_liana_bulk_%j.err

set -eo pipefail
mkdir -p scripts/figures/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/figures/figS_liana_bulk_reverse.R
