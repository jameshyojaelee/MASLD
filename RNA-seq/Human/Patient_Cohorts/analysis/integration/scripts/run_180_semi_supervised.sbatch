#!/bin/bash
#SBATCH --job-name=180_semi_supervised
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=results/novel_ml/semi_supervised/180_semi_supervised_%j.out
#SBATCH --error=results/novel_ml/semi_supervised/180_semi_supervised_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== Script 180: Semi-Supervised Learning ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

# Create output directory
mkdir -p results/novel_ml/semi_supervised

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

python scripts/180_semi_supervised.py

echo ""
echo "End: $(date)"
