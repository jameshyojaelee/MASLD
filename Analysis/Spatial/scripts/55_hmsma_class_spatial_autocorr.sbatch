#!/usr/bin/env bash
#SBATCH --job-name=squidpy
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/hmsma_class_spatial_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/hmsma_class_spatial_%j.out
set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n spatial python \
  Analysis/Spatial/scripts/55_hmsma_class_spatial_autocorr.py
echo "EXIT=$?"
