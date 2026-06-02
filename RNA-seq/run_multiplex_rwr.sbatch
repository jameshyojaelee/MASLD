#!/bin/bash
#SBATCH --job-name=multiplex_rwr
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/multiplex_rwr_%j.out
#SBATCH --error=logs/multiplex_rwr_%j.err

set -euo pipefail

echo "=== Multiplex RWR ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs

# Activate environment
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

Rscript 59_multiplex_rwr.R

echo "Done: $(date)"
