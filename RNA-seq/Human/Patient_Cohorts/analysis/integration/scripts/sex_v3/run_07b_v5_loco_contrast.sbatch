#!/bin/bash
#SBATCH --job-name=v6_07b_loco
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --array=1-7
#SBATCH --output=logs/M07b_loco_%A_%a.out
#SBATCH --error=logs/M07b_loco_%A_%a.err
# Per-fold LOCO posterior-predictive sign concordance for a given
# CONTRAST_NAME. Default array size = 7 (matches mash_vs_masl).
# disease_vs_ctrl uses 5 folds -- override --array=1-5 from CLI.
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

: "${CONTRAST_NAME:?CONTRAST_NAME must be set}"
export LOCO_FOLD="${SLURM_ARRAY_TASK_ID}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== 07b LOCO fold ${LOCO_FOLD} | CONTRAST_NAME=${CONTRAST_NAME} ==="
echo "Host: $(hostname)  date: $(date)"

Rscript 07b_v5_loco_contrast.R

echo "=== done: $(date) ==="
