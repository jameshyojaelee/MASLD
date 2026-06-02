#!/bin/bash
#SBATCH --job-name=A7_T12_comp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/A7_T12_comp_%j.out
#SBATCH --error=logs/A7_T12_comp_%j.err

set -eo pipefail
mkdir -p logs

echo "=== A7 T12: Composition decomposition ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export R_PARALLEL_SEED=42

Rscript 14.5e_T12_composition.R

echo "Done: $(date)"
