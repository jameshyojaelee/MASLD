#!/bin/bash
#SBATCH --job-name=diamond_kda
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/diamond_kda_%j.out
#SBATCH --error=logs/diamond_kda_%j.err

set -euo pipefail

mkdir -p logs

echo "=== DIAMOnD + KDA ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript RNA-seq/58_diamond_kda.R

echo ""
echo "Finished: $(date)"
