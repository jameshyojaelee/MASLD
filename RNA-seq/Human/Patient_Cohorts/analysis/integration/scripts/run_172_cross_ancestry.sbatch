#!/bin/bash
#SBATCH --job-name=172_cross_ancestry
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=results/novel_ml/cross_ancestry/172_cross_ancestry_%j.out
#SBATCH --error=results/novel_ml/cross_ancestry/172_cross_ancestry_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== Script 172: Cross-Ancestry Disease Progression Model ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

python scripts/172_cross_ancestry_progression.py

echo ""
echo "End: $(date)"
