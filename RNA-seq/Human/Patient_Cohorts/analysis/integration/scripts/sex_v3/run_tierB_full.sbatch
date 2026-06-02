#!/bin/bash
#SBATCH --job-name=v6_tierB
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --output=logs/tierB_%j.out
#SBATCH --error=logs/tierB_%j.err
# Tier B end-to-end: Modules 01 -> 02 -> 03 -> 19 for an underpowered contrast.
# Single SLURM job (sequential; modules 01-03 are short). Caller exports
# CONTRAST_NAME via sbatch --export=ALL,CONTRAST_NAME=...
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

: "${CONTRAST_NAME:?CONTRAST_NAME must be set (masl_vs_ctrl or mash_vs_ctrl)}"
case "${CONTRAST_NAME}" in
  masl_vs_ctrl|mash_vs_ctrl) ;;
  *) echo "FATAL: CONTRAST_NAME='${CONTRAST_NAME}' is not a Tier B contrast" >&2; exit 2 ;;
esac

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== Tier B end-to-end: CONTRAST_NAME=${CONTRAST_NAME} ==="
echo "Host: $(hostname)  date: $(date)"

echo "--- Module 01 (load_data) ---"
Rscript 01_load_data.R

echo "--- Module 02 (sva_factors) ---"
Rscript 02_sva_factors.R

echo "--- Module 03 (age_imputation) ---"
Rscript 03_age_imputation.R

echo "--- Module 19 (light interaction) ---"
Rscript 19_light_interaction_contrast.R

echo "=== done: $(date) ==="
