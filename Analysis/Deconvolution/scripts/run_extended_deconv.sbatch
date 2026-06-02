#!/bin/bash
#SBATCH --job-name=ext_deconv
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/deconvolution/logs/extended_deconv_%j.out

set -euo pipefail

echo "Starting extended deconvolution analysis at $(date)"

# Activate micromamba
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/deconvolution
python scripts/extended_deconv_analysis.py

echo "Completed at $(date)"
