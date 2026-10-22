#!/bin/bash
#SBATCH --job-name=r8-catalog
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/r8_catalog_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/r8_catalog_%j.err
set -euo pipefail
module load Mesa/23.1.9-GCCcore-13.2.0
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
TASK_OUTPUT="$TASK_ROOT/figures/candidates/catalog-r8-completion-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
export MASLD_PROJECT_ROOT="$TASK_ROOT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PASSPORT_MESA_MODULE=Mesa/23.1.9-GCCcore-13.2.0
cd "$TASK_ROOT"
printf '%s\n' "$TASK_OUTPUT" > "$TASK_ROOT/logs/r8_catalog_${SLURM_JOB_ID}.run_root.txt"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -B scripts/portal/complete_r8_catalog_review.py --project-root "$TASK_ROOT" --output "$TASK_OUTPUT"
