#!/bin/bash
#SBATCH --job-name=V2screenrefit
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/V2refetch/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/V2refetch/%x_%j.err

# V2 sensitivity refit of the 774-unit disease model.
#   Set VIRAL_MODE, VIRAL_PERSON_CALLS and VIRAL_OUTPUT_ROOT before submission.
set -euo pipefail
umask 027

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
: "${VIRAL_MODE:?VIRAL_MODE is required}"
: "${VIRAL_PERSON_CALLS:?VIRAL_PERSON_CALLS is required}"
: "${VIRAL_OUTPUT_ROOT:?VIRAL_OUTPUT_ROOT is required}"
test ! -e "${VIRAL_OUTPUT_ROOT}"

export MASLD_PROJECT_ROOT="${project_root}"
export VIRAL_MODE VIRAL_PERSON_CALLS VIRAL_OUTPUT_ROOT
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export LANG="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
export PATH="${rnaseq_env}/bin:${PATH}"

"${rnaseq_env}/bin/Rscript" --vanilla \
  "${project_root}/scripts/analysis/hypotheses_2026q4/V2_viral_screen/viral_positive_refit.R"
echo "VIRAL_REFIT_COMPLETE ${VIRAL_MODE} ${VIRAL_OUTPUT_ROOT}"
