#!/bin/bash
#SBATCH --job-name=hotspot_500_tcells
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/500_tcells_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/500_tcells_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

echo "[$(date)] Starting hotspot_500_tcells on $(hostname)"

micromamba run -n rapids_singlecell python \
    Analysis/SingleCell/scripts/hotspot_modules/500_prepare_tcell_subset.py

echo "[$(date)] Finished hotspot_500_tcells"
