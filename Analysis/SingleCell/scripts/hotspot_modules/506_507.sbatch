#!/bin/bash
#SBATCH --job-name=hotspot_506_507
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/506_507_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/506_507_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/hotspot_modules/506_bulk_replication.R
micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/507_olink_spatial_validation.py
