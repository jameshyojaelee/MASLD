#!/bin/bash
#SBATCH --job-name=A7_T5_age
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/A7_T5_age_%j.out
#SBATCH --error=logs/A7_T5_age_%j.err

set -eo pipefail
mkdir -p logs

echo "=== A7 T5: Sex x age confound check ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export R_PARALLEL_SEED=42

Rscript 14.5b_T5_age_confound.R

echo "Done: $(date)"
