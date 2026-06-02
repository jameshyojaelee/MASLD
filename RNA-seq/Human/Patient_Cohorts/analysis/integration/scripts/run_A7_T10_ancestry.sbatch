#!/bin/bash
#SBATCH --job-name=A7_T10_anc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/A7_T10_anc_%j.out
#SBATCH --error=logs/A7_T10_anc_%j.err

set -eo pipefail
mkdir -p logs

echo "=== A7 T10: Ancestry-proxy interaction ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export R_PARALLEL_SEED=42

Rscript 14.5c_T10_ancestry.R

echo "Done: $(date)"
