#!/bin/bash
#SBATCH --job-name=v2_smoke
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=Analysis/SingleCell/scripts/v2/logs/02_smoke_%j.log
#SBATCH --error=Analysis/SingleCell/scripts/v2/logs/02_smoke_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/scripts/v2/logs

echo "[$(date)] Job ${SLURM_JOB_ID} starting on $(hostname)"

micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/v2/02_smoke_atlas_v2.py

echo "[$(date)] DONE"
