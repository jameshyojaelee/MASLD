#!/bin/bash
#SBATCH --job-name=hotspot_smoke_chol
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/501_smoke_chol_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/501_smoke_chol_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

echo "[$(date)] Starting cholangiocyte smoke run"
micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py \
  --cell-type cholangiocytes

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/validate.py hotspot cholangiocytes

echo "[$(date)] Smoke run complete"
