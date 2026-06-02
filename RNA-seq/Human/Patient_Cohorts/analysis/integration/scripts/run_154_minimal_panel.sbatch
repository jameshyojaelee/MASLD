#!/bin/bash
#SBATCH --job-name=154_panel
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=../results/prognosis/154_panel_%j.out
#SBATCH --error=../results/prognosis/154_panel_%j.err

set -euo pipefail

echo "=== 154_minimal_prognosis_panel ==="
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"

mkdir -p ../results/prognosis

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

python 154_minimal_prognosis_panel.py

echo "Done: $(date)"
