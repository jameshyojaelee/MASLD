#!/bin/bash
#SBATCH --job-name=script70_embeddings
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/script70_%j.out
#SBATCH --error=logs/script70_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== Script 70: Classical on Embeddings (LOCO-specific embeddings) ==="
echo "Started: $(date)"
python 70_classical_on_embeddings.py 2>&1
echo "Completed: $(date)"
