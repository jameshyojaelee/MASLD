#!/usr/bin/env bash
#SBATCH --job-name=commot_validation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/commot_validation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/commot_validation_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python Analysis/Spatial/scripts/16h_commot_validation.py
