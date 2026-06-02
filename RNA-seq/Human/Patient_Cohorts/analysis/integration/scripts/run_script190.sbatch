#!/bin/bash
#SBATCH --job-name=script190_multiprogram
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/script190_%j.out
#SBATCH --error=logs/script190_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== Script 190: Multiprogram Predictor (per-fold variance ranking) ==="
echo "Started: $(date)"
python 190_multiprogram_predictor.py 2>&1
echo "Completed: $(date)"
