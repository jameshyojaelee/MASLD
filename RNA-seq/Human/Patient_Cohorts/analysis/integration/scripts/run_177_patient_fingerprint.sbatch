#!/bin/bash
#SBATCH --job-name=177_patient_fp
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=results/novel_ml/patient_fingerprint/177_patient_fingerprint_%j.out
#SBATCH --error=results/novel_ml/patient_fingerprint/177_patient_fingerprint_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== Script 177: Integrated Patient Fingerprint ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

python scripts/177_patient_fingerprint.py

echo ""
echo "End: $(date)"
