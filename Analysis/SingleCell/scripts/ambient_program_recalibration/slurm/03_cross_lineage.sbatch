#!/bin/bash
#SBATCH --job-name=amb117_lineage
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/logs/03_lineage_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/logs/03_lineage_%j.log

set -euo pipefail
umask 027

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CAND_ROOT=${MASLD_PROJECT_ROOT}/Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2
export SCRIPT_ROOT=${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/ambient_program_recalibration
export PYTHONPATH=${SCRIPT_ROOT}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONNOUSERSITE=1
python_bin=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python

echo "[lineage] host=$(hostname) job=${SLURM_JOB_ID} start=$(date -u)"
"${python_bin}" "${SCRIPT_ROOT}/03_cross_lineage.py"
echo "[lineage] done=$(date -u)"
