#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --job-name=figS09d
#SBATCH --output=logs/figS09d_%j.out
#SBATCH --error=logs/figS09d_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
${RSCRIPT} scripts/figures/figS09d_mesusie_susiex_comparison.R
