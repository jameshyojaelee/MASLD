#!/bin/bash
#SBATCH --job-name=convergence
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/46d_convergence_%j.out
#SBATCH --error=logs/46d_convergence_%j.err

set -euo pipefail
mkdir -p logs

echo "=== 46d Convergence Evidence (reframed heuristic) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

micromamba run -n rnaseq Rscript 46d_convergence_evidence.R

echo "End: $(date)"
