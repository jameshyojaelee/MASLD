#!/usr/bin/env bash
# Packed bigmem stage for MASL-vs-Healthy:
#   - 1 dream mega fit (05g) + 5 LOO folds in parallel inside one node
#   - then 1 LFC sweep
# All on one 72c/512G bigmem allocation (~1-2h end-to-end).

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/mash_loo_sweep

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -euo pipefail
cd "${SCRIPTS}"

echo "MASL stage node $(hostname)  cpu=$(nproc)  mem=$(free -g | awk '/^Mem:/ {print $2"G"}')  start $(date)"

# Dream mega + 5 LOO folds in parallel (independent, no dep between them).
# Allocation: mega=24 CPU, each LOO fold=9 CPU, total 24+5*9=69 of 72 available.
COHORTS_MH=(GSE126848 GSE130970 GSE135251 GSE162694)

( export SLURM_CPUS_PER_TASK=24
  Rscript 05g_masl_vs_healthy_dream.R \
    > "${LOGS}/stage_masl_mega.out" 2> "${LOGS}/stage_masl_mega.err"
) &
pid_mega=$!

declare -a pids_loo=()
for cohort in "${COHORTS_MH[@]}"; do
    (
        export CONTRAST="masl_vs_healthy"
        export HELD_OUT="${cohort}"
        export SLURM_CPUS_PER_TASK=9
        Rscript dream_loo_cv_contrasts.R \
            > "${LOGS}/stage_masl_loo_${cohort}.out" \
            2> "${LOGS}/stage_masl_loo_${cohort}.err"
    ) &
    pids_loo+=($!)
done

set +e
wait ${pid_mega}; rc_mega=$?
rc_loo=0
for p in "${pids_loo[@]}"; do wait "${p}" || rc_loo=$?; done
set -e

echo "MASL stage mega rc=${rc_mega}  loo rc=${rc_loo}  $(date)"
if [[ ${rc_mega} -ne 0 || ${rc_loo} -ne 0 ]]; then
    echo "MASL stage FAILED before sweep"
    exit 1
fi

echo ""
echo "MASL sweep starting $(date)"
( export CONTRAST="masl_vs_healthy" SLURM_CPUS_PER_TASK=4
  Rscript lfc_sweep_loo.R \
    > "${LOGS}/stage_masl_sweep.out" \
    2> "${LOGS}/stage_masl_sweep.err"
)
echo "MASL stage OK $(date)"
