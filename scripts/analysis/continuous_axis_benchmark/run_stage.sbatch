#!/bin/bash
# Generic runner for any stage of CONTINUOUS-AXIS-BENCHMARK-v1.
#
# Required: CAB_OUT_ROOT (an existing root created by 01), CAB_SCRIPT (basename).
# Optional: CAB_ARM, CAB_MATRIX, CAB_N_SPLIT -- passed through to the script.
#
# Threads are pinned to 1 deliberately. Threaded BLAS buys nothing for the row-wise
# work here, and preprocessCore already demonstrated that libraries probing the
# node core count rather than the cgroup allocation abort under SLURM.
#SBATCH --job-name=cabstage
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/cab_%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/cab_%x_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"

: "${CAB_OUT_ROOT:?Set CAB_OUT_ROOT to the output root created by 01}"
: "${CAB_SCRIPT:?Set CAB_SCRIPT to a script basename in scripts/analysis/continuous_axis_benchmark}"
test -d "${CAB_OUT_ROOT}" || { echo "CAB_OUT_ROOT does not exist: ${CAB_OUT_ROOT}" >&2; exit 1; }
test -f "${CAB_OUT_ROOT}/EXEMPTION.json" || { echo "No EXEMPTION.json in the output root" >&2; exit 1; }

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

mkdir -p "${PROJECT_ROOT}/logs" "${CAB_OUT_ROOT}/logs"
cd "${PROJECT_ROOT}"

STEM="$(basename "${CAB_SCRIPT}" .R)"
TAG="${STEM}${CAB_ARM:+_${CAB_ARM}}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "scripts/analysis/continuous_axis_benchmark/${CAB_SCRIPT}" \
  2>&1 | tee "${CAB_OUT_ROOT}/logs/${TAG}.log"

printf '%s\n' "${CAB_OUT_ROOT}"
