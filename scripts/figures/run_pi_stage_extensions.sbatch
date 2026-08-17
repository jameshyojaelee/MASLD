#!/bin/bash
#SBATCH --job-name=pi_stage_ext
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/pi_stage_ext_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/pi_stage_ext_%j.err

set -euo pipefail
umask 027

: "${FIGURE_CANDIDATE_ROOT:?Set a new timestamped candidate path}"
: "${STAGE_RELEASE_ROOT:?Set the validated fragment-native stage release}"
test ! -e "${FIGURE_CANDIDATE_ROOT}"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PATH="${RNASEQ_ENV}/bin:${PATH}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "${PROJECT_ROOT}/scripts/figures/build_pi_stage_extensions.R"
