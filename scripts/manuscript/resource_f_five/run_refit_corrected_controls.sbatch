#!/bin/bash
#SBATCH --job-name=fcorr_refit
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/resource_f_five/fcorr_%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/resource_f_five/fcorr_%x_%j.err

# Corrected GSE130970 control refit of the five-cohort F_five model.
#   REFIT_MODE=fit        REFIT_OUTPUT_ROOT=.../bulk-corrected-controls-<UTC>
#   REFIT_MODE=reproduce  REFIT_OUTPUT_ROOT=.../bulk-corrected-controls-<UTC>-reproduction
#                         REFIT_SOURCE_ROOT=.../bulk-corrected-controls-<UTC>
# Both roots must not exist; the driver refuses to write into an existing root.
set -euo pipefail
umask 027

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
: "${REFIT_OUTPUT_ROOT:?REFIT_OUTPUT_ROOT is required}"
export REFIT_MODE="${REFIT_MODE:-fit}"
test ! -e "${REFIT_OUTPUT_ROOT}"

export MASLD_PROJECT_ROOT="${project_root}"
export REFIT_OUTPUT_ROOT
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export LANG="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
export PATH="${rnaseq_env}/bin:${PATH}"

"${rnaseq_env}/bin/Rscript" --vanilla \
  "${project_root}/scripts/manuscript/resource_f_five/refit_corrected_controls.R"
echo "FCORR_${REFIT_MODE^^}_COMPLETE ${REFIT_OUTPUT_ROOT}"
