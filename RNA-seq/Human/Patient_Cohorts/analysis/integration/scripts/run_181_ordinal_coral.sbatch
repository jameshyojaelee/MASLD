#!/bin/bash
#SBATCH --job-name=181_ordinal_coral
#SBATCH --output=scripts/logs/181_ordinal_coral_%j.out
#SBATCH --error=scripts/logs/181_ordinal_coral_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00

set -euo pipefail

echo "[181] Starting ordinal CORAL v2 — $(date)"
echo "[181] Job ID: ${SLURM_JOB_ID}"
echo "[181] Node: $(hostname)"

# Activate environment (torch available in spatial)
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

python scripts/181_ordinal_coral_v2.py

echo "[181] Completed — $(date)"
