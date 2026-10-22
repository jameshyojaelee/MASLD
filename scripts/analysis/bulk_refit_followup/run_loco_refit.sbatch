#!/bin/bash
#SBATCH --job-name=loco_refit
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/bulk_refit_followup/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/bulk_refit_followup/%x_%j.err

# Leave-one-cohort-out refit of the 774-unit primary_source_controls model.
#   sbatch --job-name=loco_<COHORT> \
#     --export=ALL,LOCO_COHORT=<COHORT|NONE>,LOCO_OUTPUT_ROOT=<new candidate dir> run_loco_refit.sbatch
set -euo pipefail
umask 027

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
: "${LOCO_COHORT:?LOCO_COHORT is required}"
: "${LOCO_OUTPUT_ROOT:?LOCO_OUTPUT_ROOT is required}"
test ! -e "${LOCO_OUTPUT_ROOT}"

export MASLD_PROJECT_ROOT="${project_root}"
export LOCO_COHORT LOCO_OUTPUT_ROOT
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export LANG="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
export PATH="${rnaseq_env}/bin:${PATH}"

"${rnaseq_env}/bin/Rscript" --vanilla \
  "${project_root}/scripts/analysis/bulk_refit_followup/loco_refit.R"
echo "LOCO_COMPLETE ${LOCO_COHORT} ${LOCO_OUTPUT_ROOT}"
