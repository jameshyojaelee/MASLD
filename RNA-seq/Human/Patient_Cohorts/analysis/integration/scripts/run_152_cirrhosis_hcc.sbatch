#!/bin/bash
#SBATCH --job-name=152_cirrhosis_hcc
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

echo "=== Script 152: Cirrhosis + HCC Molecular Risk ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

python scripts/152_cirrhosis_hcc_risk.py

echo ""
echo "End: $(date)"
