#!/usr/bin/env bash
# Stage 1 (cpu, PRJNA-excluded): 4 dream mega fits in parallel inside one
# cpu node. Each fit uses ~7 CPU with BLAS pinned (4 procs × 7 = 28 CPU).

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

echo "Stage 1 all-megas (PRJNA excl) node $(hostname)  cpu=$(nproc)  mem=$(free -g | awk '/^Mem:/ {print $2"G"}')  start $(date)"

( export SLURM_CPUS_PER_TASK=7 MASH_DEF=borderline_grouped
  Rscript 05e_mash_vs_masl_dream_strict.R \
    > "${LOGS}/stage1_mm_primary.out" 2> "${LOGS}/stage1_mm_primary.err"
) &
pid_a=$!

( export SLURM_CPUS_PER_TASK=7 MASH_DEF=strict
  Rscript 05e_mash_vs_masl_dream_strict.R \
    > "${LOGS}/stage1_mm_strict.out" 2> "${LOGS}/stage1_mm_strict.err"
) &
pid_b=$!

( export SLURM_CPUS_PER_TASK=7 MASH_DEF=borderline_grouped
  Rscript 05f_mash_vs_healthy_dream.R \
    > "${LOGS}/stage1_mh_primary.out" 2> "${LOGS}/stage1_mh_primary.err"
) &
pid_c=$!

( export SLURM_CPUS_PER_TASK=7 MASH_DEF=strict
  Rscript 05f_mash_vs_healthy_dream.R \
    > "${LOGS}/stage1_mh_strict.out" 2> "${LOGS}/stage1_mh_strict.err"
) &
pid_d=$!

set +e
wait ${pid_a}; rc_a=$?
wait ${pid_b}; rc_b=$?
wait ${pid_c}; rc_c=$?
wait ${pid_d}; rc_d=$?
set -e
echo "Stage 1 all rc: mm_primary=${rc_a} mm_strict=${rc_b} mh_primary=${rc_c} mh_strict=${rc_d}  $(date)"
if [[ ${rc_a} -ne 0 || ${rc_b} -ne 0 || ${rc_c} -ne 0 || ${rc_d} -ne 0 ]]; then
    echo "Stage 1 FAILED"; exit 1
fi
echo "Stage 1 all OK $(date)"
