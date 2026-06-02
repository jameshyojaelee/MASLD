#!/bin/bash
#SBATCH --job-name=scdrs_mg
#SBATCH --partition=cpu
#SBATCH --mem=160G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --array=0-18
#SBATCH --output=scripts/figures/logs/scdrs_mg_%A_%a.out
#SBATCH --error=scripts/figures/logs/scdrs_mg_%A_%a.err

# Score one MAGMA-anchored trait on the 25% subsample atlas per array task.
# 19 GWAS, %10 concurrent under nslab QoS on cpu partition.

set -eo pipefail
mkdir -p scripts/figures/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

export SCDRS_TRAIT_INDEX=$SLURM_ARRAY_TASK_ID
python scripts/figures/run_scdrs_magma_25pct.py
