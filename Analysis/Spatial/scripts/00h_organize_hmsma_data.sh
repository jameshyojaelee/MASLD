#!/bin/bash
#SBATCH --job-name=organize_hmsma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=../logs/organize_hmsma_%j.out
#SBATCH --error=../logs/organize_hmsma_%j.err

# 00h: Organize HMSMA processed data into SpaceRanger-compatible structure
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd "$(dirname "$0")"
python 00h_organize_hmsma_data.py "$@"
