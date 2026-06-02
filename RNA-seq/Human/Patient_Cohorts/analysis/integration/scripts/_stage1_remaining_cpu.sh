#!/usr/bin/env bash
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/mash_loo_sweep

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

set -euo pipefail
cd "${SCRIPTS}"

echo "Stage 1 remaining (cpu) node $(hostname)  cpu=$(nproc)  mem=$(free -g | awk '/^Mem:/ {print $2"G"}')  start $(date)"

( export SLURM_CPUS_PER_TASK=12 MASH_DEF=borderline_grouped
  Rscript 05f_mash_vs_healthy_dream.R \
    > "${LOGS}/stage1cpu_mh_primary.out" 2> "${LOGS}/stage1cpu_mh_primary.err"
) &
pid_a=$!

( export SLURM_CPUS_PER_TASK=12
  Rscript 05e_mash_vs_masl_dream_strict.R \
    > "${LOGS}/stage1cpu_mm_strict.out" 2> "${LOGS}/stage1cpu_mm_strict.err"
) &
pid_b=$!

set +e
wait ${pid_a}; rc_a=$?
wait ${pid_b}; rc_b=$?
set -e
echo "Stage 1 remaining rc: mh_primary=${rc_a} mm_strict=${rc_b}  $(date)"
[[ ${rc_a} -eq 0 && ${rc_b} -eq 0 ]] || { echo "FAILED"; exit 1; }
echo "Stage 1 remaining OK $(date)"
