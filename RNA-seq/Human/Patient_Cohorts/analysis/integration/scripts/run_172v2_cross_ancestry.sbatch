#!/bin/bash
#SBATCH --job-name=172v2_cross_ancestry
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/172v2_cross_ancestry_%j.out
#SBATCH --error=logs/172v2_cross_ancestry_%j.err

set -euo pipefail

echo "=== Script 172v2: Cross-Ancestry Staging Classifier (Fixed) ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "CPUs:   ${SLURM_CPUS_PER_TASK}"
echo "Start:  $(date)"
echo ""

# Activate spatial env (avoids scipy/numpy ABI conflict in rnaseq env)
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

python scripts/172v2_cross_ancestry_fixed.py

echo ""
echo "=== Done: $(date) ==="
