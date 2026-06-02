#!/bin/bash
#SBATCH --job-name=v6_19_lightI
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --output=logs/19_%j.out
#SBATCH --error=logs/19_%j.err
# Tier B — light fixed-effect sex × disease interaction for an underpowered
# contrast (masl_vs_ctrl or mash_vs_ctrl). One dream fit + per-cohort Cochran's
# Q heterogeneity. ~30 min wall-clock on cpu (no permutation / bootstrap /
# calibration). Caller must export CONTRAST_NAME via sbatch --export.
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
echo "=== Tier B run: CONTRAST_NAME=${CONTRAST_NAME} ==="
echo "Host: $(hostname)  date: $(date)"
Rscript 19_light_interaction_contrast.R
echo "=== done: $(date) ==="
