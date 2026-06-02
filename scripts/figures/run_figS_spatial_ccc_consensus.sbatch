#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=2 --mem=10G --time=48:00:00
#SBATCH --job-name=figS_I1 --output=scripts/figures/logs/figS_I1_%j.out --error=scripts/figures/logs/figS_I1_%j.err
set -eo pipefail; mkdir -p scripts/figures/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript scripts/figures/figS_spatial_ccc_consensus.R
