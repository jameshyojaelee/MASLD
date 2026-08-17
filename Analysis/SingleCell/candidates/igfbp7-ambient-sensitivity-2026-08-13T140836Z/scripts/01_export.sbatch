#!/bin/bash
#SBATCH --job-name=ambient
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/01_export_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/01_export_%j.log

set -euo pipefail
umask 027

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CAND_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONNOUSERSITE=1

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID} start=$(date -u)"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python \
  "${CAND_ROOT}/scripts/01_export_counts.py"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -m pip freeze \
  > "${CAND_ROOT}/results/01_pip_freeze.txt"
echo "[sbatch] DONE $(date -u)"
