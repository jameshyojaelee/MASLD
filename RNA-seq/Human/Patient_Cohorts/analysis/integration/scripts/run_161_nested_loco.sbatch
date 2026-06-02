#!/bin/bash
#SBATCH --job-name=161_nested_loco
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=results/prognosis_v2/161_nested_loco_%j.out
#SBATCH --error=results/prognosis_v2/161_nested_loco_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== Script 161: Nested LOCO Classifier ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

python scripts/161_nested_loco_classifier.py

echo ""
echo "End: $(date)"
