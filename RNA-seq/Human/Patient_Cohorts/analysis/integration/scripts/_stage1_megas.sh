#!/usr/bin/env bash
# Stage 1: 3 dream mega fits in parallel inside one bigmem node.
# Invoked by run_mash_loo_sweep.sh; not meant to be run directly.

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/mash_loo_sweep

# micromamba activate scripts reference unbound vars under `set -u`; defer
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -euo pipefail
cd "${SCRIPTS}"

echo "Stage 1 node $(hostname)  cpu=$(nproc)  mem=$(free -g | awk '/^Mem:/ {print $2"G"}')  start $(date)"

( export SLURM_CPUS_PER_TASK=24 MASH_DEF=borderline_grouped
  Rscript 05f_mash_vs_healthy_dream.R \
    > "${LOGS}/stage1_mh_primary.out" 2> "${LOGS}/stage1_mh_primary.err"
) &
pid_a=$!

( export SLURM_CPUS_PER_TASK=24 MASH_DEF=strict
  Rscript 05f_mash_vs_healthy_dream.R \
    > "${LOGS}/stage1_mh_strict.out" 2> "${LOGS}/stage1_mh_strict.err"
) &
pid_b=$!

( export SLURM_CPUS_PER_TASK=24
  Rscript 05e_mash_vs_masl_dream_strict.R \
    > "${LOGS}/stage1_mm_strict.out" 2> "${LOGS}/stage1_mm_strict.err"
) &
pid_c=$!

set +e
wait ${pid_a}; rc_a=$?
wait ${pid_b}; rc_b=$?
wait ${pid_c}; rc_c=$?
set -e

echo "Stage 1 exit codes: mh_primary=${rc_a} mh_strict=${rc_b} mm_strict=${rc_c}"
if [[ ${rc_a} -ne 0 || ${rc_b} -ne 0 || ${rc_c} -ne 0 ]]; then
    echo "Stage 1 FAILED"
    exit 1
fi
echo "Stage 1 OK $(date)"
