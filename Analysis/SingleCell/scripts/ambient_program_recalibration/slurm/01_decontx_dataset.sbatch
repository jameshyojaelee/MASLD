#!/bin/bash
#SBATCH --job-name=amb117_decontx
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/logs/01_decontx_%x_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/logs/01_decontx_%x_%j.log

set -euo pipefail
umask 027

: "${DECONTX_ONLY_DATASET:?DECONTX_ONLY_DATASET is required}"
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CAND_ROOT=${MASLD_PROJECT_ROOT}/Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2
export SCRIPT_ROOT=${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/ambient_program_recalibration
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "[decontX] host=$(hostname) job=${SLURM_JOB_ID} dataset=${DECONTX_ONLY_DATASET} start=$(date -u)"
micromamba run -n decontx Rscript "${SCRIPT_ROOT}/01_decontx_all117.R"
echo "[decontX] done=$(date -u)"
