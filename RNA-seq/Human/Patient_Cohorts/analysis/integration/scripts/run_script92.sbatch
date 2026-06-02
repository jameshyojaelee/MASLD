#!/bin/bash
#SBATCH --job-name=script92_multitask
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/script92_%j.out
#SBATCH --error=logs/script92_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

echo "=== Script 92: Multi-task Semi-supervised (per-fold genes) ==="
echo "Started: $(date)"
python 92_multitask_semisupervised.py 2>&1
echo "Completed: $(date)"
