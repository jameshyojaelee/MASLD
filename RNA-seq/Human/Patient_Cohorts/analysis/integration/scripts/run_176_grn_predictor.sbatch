#!/bin/bash
#SBATCH --job-name=176_grn_predictor
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=results/novel_ml/grn_predictor/176_grn_predictor_%j.out
#SBATCH --error=results/novel_ml/grn_predictor/176_grn_predictor_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== Script 176: Cell-Type-Resolved GRN Predictor ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

python scripts/176_grn_predictor.py

echo ""
echo "End: $(date)"
