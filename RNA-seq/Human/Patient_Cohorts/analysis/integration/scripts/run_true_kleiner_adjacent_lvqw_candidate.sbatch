#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=logs/limma_true_kleiner_%j.out
#SBATCH --error=logs/limma_true_kleiner_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
OUTPUT_ROOT="${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates/program-context-v2-candidate-2026-08-07/fibrosis-adjacent-true-kleiner-lvqw-v1"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"

# Keep threaded BLAS/OpenMP libraries within the requested SLURM allocation.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

if [[ -e "${OUTPUT_ROOT}" || -L "${OUTPUT_ROOT}" ]]; then
  echo "Refusing existing candidate output: ${OUTPUT_ROOT}" >&2
  exit 1
fi

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
"${RNASEQ_ENV}/bin/Rscript" "${SCRIPT_DIR}/14h_true_kleiner_adjacent_lvqw_candidate.R"
"${RNASEQ_ENV}/bin/python3" "${SCRIPT_DIR}/validate_true_kleiner_adjacent_lvqw_candidate.py" \
  --project-root "${PROJECT_ROOT}" \
  --output-root "${OUTPUT_ROOT}"

test -f "${OUTPUT_ROOT}/READY"
echo "READY: ${OUTPUT_ROOT}/READY"
