#!/usr/bin/env bash
#SBATCH --job-name=multisp_validation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/multisp_validation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/multisp_validation_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial_multiomics

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python Analysis/Spatial/scripts/13d_multisp_validation.py
