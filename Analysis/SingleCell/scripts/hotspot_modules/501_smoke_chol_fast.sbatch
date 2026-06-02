#!/bin/bash
#SBATCH --job-name=hotspot_smoke_fast
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/501_smoke_chol_fast_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/501_smoke_chol_fast_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

export HOTSPOT_MAX_GENES=2000
export HOTSPOT_JOBS=16

# Write to a separate output dir so it doesn't collide with the slow smoke
export HOTSPOT_OUTDIR_SUFFIX="_fast"

echo "[$(date)] Starting tuned cholangiocyte smoke (max_genes=2000, jobs=16)"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py \
  --cell-type cholangiocytes

echo "[$(date)] Tuned smoke complete"
