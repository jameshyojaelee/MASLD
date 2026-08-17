#!/bin/bash
#SBATCH --job-name=cab01
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/cab01_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/cab01_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"

: "${CAB_OUT_ROOT:?Set CAB_OUT_ROOT to a new timestamped output root}"
test ! -e "${CAB_OUT_ROOT}"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONNOUSERSITE=1
# Pin to one thread. Threaded BLAS gives no speedup for the row-wise work here and
# libraries that probe the node core count rather than the cgroup allocation abort
# under SLURM (preprocessCore did exactly that with pthread_create EINVAL).
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

mkdir -p "${PROJECT_ROOT}/logs" "${CAB_OUT_ROOT}/logs"
cd "${PROJECT_ROOT}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  scripts/analysis/continuous_axis_benchmark/01_build_substrate.R \
  2>&1 | tee "${CAB_OUT_ROOT}/logs/01_substrate.log"

cp -a scripts/analysis/continuous_axis_benchmark/*.R "${CAB_OUT_ROOT}/logs/"
printf '%s\n' "${CAB_OUT_ROOT}"
