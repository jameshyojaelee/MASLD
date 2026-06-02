#!/bin/bash
#SBATCH --job-name=v6_tierB_bun
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --output=logs/tierB_bun_%j.out
#SBATCH --error=logs/tierB_bun_%j.err
# Tier B bundle: run BOTH masl_vs_ctrl and mash_vs_ctrl pipelines in parallel
# inside one 512G allocation (each pipeline peaks ~180G; 2x180=360G fits in
# 512G with 152G buffer). Each gets 8 CPUs of the 16 available.
set -e
SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

run_one() {
  local CN=$1
  local LOG="${SCRIPT_DIR}/logs/tierB_bun_${SLURM_JOB_ID}_${CN}.log"
  {
    echo "=== Tier B inner-parallel: ${CN}  host=$(hostname)  $(date) ==="
    CONTRAST_NAME=${CN} OMP_NUM_THREADS=8 SLURM_CPUS_PER_TASK=8 \
      Rscript 01_load_data.R                && \
    CONTRAST_NAME=${CN} OMP_NUM_THREADS=8 SLURM_CPUS_PER_TASK=8 \
      Rscript 02_sva_factors.R              && \
    CONTRAST_NAME=${CN} OMP_NUM_THREADS=8 SLURM_CPUS_PER_TASK=8 \
      Rscript 03_age_imputation.R           && \
    CONTRAST_NAME=${CN} OMP_NUM_THREADS=8 SLURM_CPUS_PER_TASK=8 \
      Rscript 19_light_interaction_contrast.R
    echo "=== Tier B done: ${CN}  $(date) ==="
  } >> "$LOG" 2>&1
}

export -f run_one
export SCRIPT_DIR

echo "=== Starting both Tier B pipelines in parallel ==="
date
run_one masl_vs_ctrl &
PID1=$!
run_one mash_vs_ctrl &
PID2=$!

wait $PID1
RC1=$?
wait $PID2
RC2=$?

echo "=== Tier B bundle done.  masl_vs_ctrl rc=$RC1  mash_vs_ctrl rc=$RC2 ==="
date
[[ $RC1 -eq 0 && $RC2 -eq 0 ]] || exit 1
