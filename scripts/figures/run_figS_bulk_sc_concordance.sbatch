#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --job-name=figS_bsc_conc
#SBATCH --output=scripts/figures/logs/figS_bulk_sc_concordance_%j.out
#SBATCH --error=scripts/figures/logs/figS_bulk_sc_concordance_%j.err

set -eo pipefail
mkdir -p scripts/figures/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/figures/figS_bulk_sc_concordance.R
