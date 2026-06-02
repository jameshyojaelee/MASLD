#!/bin/bash
#SBATCH --job-name=150_prognosis_feat
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=../results/prognosis/150_prognosis_%j.out
#SBATCH --error=../results/prognosis/150_prognosis_%j.err

set -euo pipefail

echo "=== 150_prognosis_feature_engineering ==="
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"

# Create output directory for logs
mkdir -p ../results/prognosis

# Activate environment
eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

python 150_prognosis_feature_engineering.py

echo "Done: $(date)"
