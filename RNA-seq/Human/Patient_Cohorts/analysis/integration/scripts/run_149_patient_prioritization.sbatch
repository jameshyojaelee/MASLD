#!/bin/bash
#SBATCH --job-name=149_patient_targets
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/149_patient_target_prioritization_%j.log
#SBATCH --error=logs/149_patient_target_prioritization_%j.log

set -euo pipefail

echo "=== Script 149: Per-Patient Therapeutic Target Prioritization ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration
mkdir -p logs

python scripts/149_patient_target_prioritization.py

echo ""
echo "End: $(date)"
