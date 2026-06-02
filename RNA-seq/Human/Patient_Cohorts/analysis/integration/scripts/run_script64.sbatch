#!/bin/bash
#SBATCH --job-name=script64_model_sweep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/script64_%j.out
#SBATCH --error=logs/script64_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== Script 64: Plan1 Model Sweep (per-fold genes) ==="
echo "Started: $(date)"
python 64_plan1_model_sweep.py 2>&1
echo "Completed: $(date)"
