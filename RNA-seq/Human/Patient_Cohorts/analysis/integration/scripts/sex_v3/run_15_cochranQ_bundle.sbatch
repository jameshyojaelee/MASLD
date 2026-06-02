#!/bin/bash
#SBATCH --job-name=v6_P6_bun
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=12:00:00
#SBATCH --output=logs/P6_bun_%j.out
#SBATCH --error=logs/P6_bun_%j.err
# P6 cochranQ bundle: run all 5 fold-out cohort tasks in parallel inside one
# 512G allocation. Each task peaks ~64G; 5x64=320G fits in 512G.
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
N_TASKS="${N_TASKS:-5}"

: "${CONTRAST_NAME:?CONTRAST_NAME must be set (passed via sbatch --export=ALL,CONTRAST_NAME=...)}"

run_task() {
  local k=$1
  local LOG="${SCRIPT_DIR}/logs/P6_bun_${SLURM_JOB_ID}_t${k}.log"
  {
    echo "=== P6 task ${k}  contrast=${CONTRAST_NAME}  host=$(hostname)  $(date) ==="
    CONTRAST_NAME=${CONTRAST_NAME} SLURM_ARRAY_TASK_ID=${k} \
      OMP_NUM_THREADS=3 SLURM_CPUS_PER_TASK=3 \
      Rscript 15_cochranQ_per_cohort.R
    echo "=== P6 task ${k} done  $(date) ==="
  } >> "$LOG" 2>&1
}
export CONTRAST_NAME

export -f run_task
export SCRIPT_DIR

echo "=== P6 cochranQ bundle: ${N_TASKS} tasks in parallel ==="
date
PIDS=()
for k in $(seq 1 ${N_TASKS}); do
  run_task $k &
  PIDS+=($!)
done
RC=0
for pid in "${PIDS[@]}"; do
  wait $pid || RC=$?
done
echo "=== P6 bundle done.  combined rc=$RC ==="
date
exit $RC
