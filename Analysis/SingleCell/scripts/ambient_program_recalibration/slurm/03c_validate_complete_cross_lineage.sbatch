#!/bin/bash
#SBATCH --job-name=amb117_lineage6val
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2/logs/03c_validate_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2/logs/03c_validate_%j.log

set -euo pipefail
umask 027

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CROSS_CAND_ROOT=${MASLD_PROJECT_ROOT}/Analysis/SingleCell/candidates/cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2
export SCRIPT_ROOT=${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/ambient_program_recalibration
export CAND_ROOT=${MASLD_PROJECT_ROOT}/Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2
export PYTHONPATH=${SCRIPT_ROOT}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONNOUSERSITE=1
python_bin=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python

echo "[lineage6-validator] host=$(hostname) job=${SLURM_JOB_ID} start=$(date -u)"
"${python_bin}" "${SCRIPT_ROOT}/03c_validate_complete_cross_lineage.py"
echo "[lineage6-validator] done=$(date -u)"
