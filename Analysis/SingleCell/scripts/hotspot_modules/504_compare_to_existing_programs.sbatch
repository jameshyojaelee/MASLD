#!/bin/bash
#SBATCH --job-name=hotspot_504
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/504_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/504_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs
micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/504_compare_to_existing_programs.py
