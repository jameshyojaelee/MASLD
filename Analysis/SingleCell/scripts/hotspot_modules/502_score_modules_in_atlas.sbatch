#!/bin/bash
#SBATCH --job-name=hotspot_502
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/502_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/502_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs
micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/502_score_modules_in_atlas.py
