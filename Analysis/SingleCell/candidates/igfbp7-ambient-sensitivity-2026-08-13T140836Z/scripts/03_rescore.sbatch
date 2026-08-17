#!/bin/bash
#SBATCH --job-name=ambient
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/03_rescore_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/03_rescore_%j.log

set -euo pipefail
umask 027

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CAND_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID} start=$(date -u)"
"${PY}" "${CAND_ROOT}/scripts/03_rescore.py"
"${PY}" "${CAND_ROOT}/scripts/04_diagnostics.py"
"${PY}" -m pip freeze > "${CAND_ROOT}/results/03_pip_freeze.txt"
echo "[sbatch] DONE $(date -u)"
