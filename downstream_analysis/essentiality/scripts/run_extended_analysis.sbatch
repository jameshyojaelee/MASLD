#!/bin/bash
#SBATCH --job-name=extended_essentiality
#SBATCH --output=logs/extended_%j.out
#SBATCH --error=logs/extended_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu

mkdir -p logs

# Activate Environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/.mamba/music_deconv

# Environment fixes
export LD_LIBRARY_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/.mamba/music_deconv/lib:$LD_LIBRARY_PATH
export MPLBACKEND=Agg
export PYTHONNOUSERSITE=1

echo "Starting extended analysis..."
python3 extended_essentiality_analysis.py
echo "Done."
