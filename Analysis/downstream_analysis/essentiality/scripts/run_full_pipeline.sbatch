#!/bin/bash
#SBATCH --job-name=full_essentiality
#SBATCH --output=logs/full_essentiality_%j.out
#SBATCH --error=logs/full_essentiality_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu

mkdir -p logs

# Activate Environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/.mamba/music_deconv

# Environment fixes
export LD_LIBRARY_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/RNA-seq/.mamba/music_deconv/lib:$LD_LIBRARY_PATH
export MPLBACKEND=Agg
export PYTHONNOUSERSITE=1

echo "Starting FULL essentiality pipeline..."
python3 run_essentiality_pipeline.py
echo "Done."
