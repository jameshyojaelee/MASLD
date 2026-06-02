#!/bin/bash
#SBATCH --job-name=magma_gs
#SBATCH --partition=cpu
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --output=scripts/figures/logs/magma_gs_%j.out
#SBATCH --error=scripts/figures/logs/magma_gs_%j.err

# Convert per-GWAS MAGMA .genes.out → unified MAGMA_anchored_scdrs.gs.
# Runs after the magma_ga array completes via --dependency=afterok.

set -eo pipefail
mkdir -p scripts/figures/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

python scripts/figures/magma_to_gs.py
