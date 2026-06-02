#!/bin/bash
#SBATCH --job-name=v6_loco_bun
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=24:00:00
#SBATCH --output=logs/loco_bun_%j.out
#SBATCH --error=logs/loco_bun_%j.err
# LOCO bundle: run multiple folds in parallel inside one 512G allocation.
# Caller passes FOLD_START + FOLD_END (1..7 for mash_vs_masl, 1..5 for
# disease_vs_ctrl). Each fold gets 4 CPUs of the 16 available.
# 3 folds in parallel x 150G peak = 450G; fits in 512G.
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

: "${CONTRAST_NAME:?CONTRAST_NAME must be set}"
: "${FOLD_START:?FOLD_START must be set (1-N)}"
: "${FOLD_END:?FOLD_END must be set (1-N)}"

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

run_fold() {
  local k=$1
  local LOG="${SCRIPT_DIR}/logs/loco_bun_${SLURM_JOB_ID}_fold${k}.log"
  {
    echo "=== fold ${k}  ${CONTRAST_NAME}  host=$(hostname)  $(date) ==="
    LOCO_FOLD=${k} OMP_NUM_THREADS=4 SLURM_CPUS_PER_TASK=4 \
      Rscript 07b_v5_loco_contrast.R
    echo "=== fold ${k} done  $(date) ==="
  } >> "$LOG" 2>&1
}

export -f run_fold
export SCRIPT_DIR CONTRAST_NAME

echo "=== LOCO bundle: folds ${FOLD_START}..${FOLD_END}  CONTRAST=${CONTRAST_NAME} ==="
date
PIDS=()
for k in $(seq ${FOLD_START} ${FOLD_END}); do
  run_fold $k &
  PIDS+=($!)
done
RC=0
for pid in "${PIDS[@]}"; do
  wait $pid || RC=$?
done
echo "=== LOCO bundle done.  combined rc=$RC ==="
date
exit $RC
