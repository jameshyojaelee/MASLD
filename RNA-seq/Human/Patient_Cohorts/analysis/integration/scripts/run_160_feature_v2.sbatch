#!/bin/bash
#SBATCH --job-name=160_feat_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=../results/prognosis_v2/160_feat_v2_%j.out
#SBATCH --error=../results/prognosis_v2/160_feat_v2_%j.err

set -euo pipefail

echo "=== 160_feature_engineering_v2 ==="
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"

# Create output directory for logs
mkdir -p ../results/prognosis_v2

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

python 160_feature_engineering_v2.py

echo "Done: $(date)"
