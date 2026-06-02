#!/bin/bash
#SBATCH --job-name=plasma_225ab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/plasma_225ab_%j.out
#SBATCH --error=logs/plasma_225ab_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== Script 225a: Plasma Sweep Linear (per-fold imputation) ==="
echo "Started: $(date)"
python 225a_plasma_sweep_linear.py 2>&1
echo "Completed 225a: $(date)"

echo "=== Script 225b: Plasma Sweep Boosting (per-fold imputation) ==="
echo "Started: $(date)"
python 225b_plasma_sweep_boosting.py 2>&1
echo "Completed 225b: $(date)"

echo "All done: $(date)"
