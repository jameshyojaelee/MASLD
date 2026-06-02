#!/bin/bash
#SBATCH --job-name=essentiality_analysis
#SBATCH --output=logs/essentiality_%j.out
#SBATCH --error=logs/essentiality_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

# Ensure log directory exists
mkdir -p logs

# Activate Environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/.mamba/music_deconv

# Fix libstdc++ issue by pointing to the conda env's lib directory
export LD_LIBRARY_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/.mamba/music_deconv/lib:$LD_LIBRARY_PATH

# Run Script
python3 RNA-seq/essentiality/scripts/analyze_essentiality.py
