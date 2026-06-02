#!/usr/bin/env bash
#SBATCH --job-name=atac_pseudobulk_proj
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/atac_pseudobulk_proj_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/atac_pseudobulk_proj_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial_multiomics

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
python Analysis/Spatial/scripts/13a_atac_pseudobulk_projection.py
