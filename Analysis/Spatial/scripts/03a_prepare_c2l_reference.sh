#!/usr/bin/env bash
#SBATCH --job-name=c2l_reference
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/c2l_reference_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/c2l_reference_%j.err
set -euo pipefail

module load CUDA/12.1.1

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python Analysis/Spatial/scripts/03a_prepare_c2l_reference.py
