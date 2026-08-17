#!/bin/bash
#SBATCH --job-name=ambient
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/06_assemble_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/06_assemble_%j.log

set -euo pipefail
umask 027
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CAND_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z
export PYTHONNOUSERSITE=1
echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID} start=$(date -u)"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python \
  "${CAND_ROOT}/scripts/06_assemble.py"
echo "[sbatch] DONE $(date -u)"
