#!/bin/bash
#SBATCH --job-name=plasma_225cg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/plasma_225cg_%j.out
#SBATCH --error=logs/plasma_225cg_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== Script 225c: Plasma Sweep Neural (per-fold imputation + per-fold TabPFN features) ==="
echo "Started: $(date)"
python 225c_plasma_sweep_neural.py 2>&1
echo "Completed 225c: $(date)"

echo "=== Script 225g: Plasma LOEO Baselines (per-fold imputation) ==="
echo "Started: $(date)"
python 225g_plasma_loeo_baselines.py 2>&1
echo "Completed 225g: $(date)"

echo "All done: $(date)"
