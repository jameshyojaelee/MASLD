#!/usr/bin/env bash
# MASL-vs-Healthy packed cpu version: mega + 5 LOO folds in parallel, then sweep.
# Fits in 28 CPU / 210G on a standard cpu node.

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

echo "MASL cpu stage  node $(hostname)  cpu=$(nproc)  mem=$(free -g | awk '/^Mem:/ {print $2"G"}')  start $(date)"

COHORTS_MH=(GSE126848 GSE130970 GSE135251 GSE162694)

# Mega (6 CPU) + 5 LOO folds × 4 CPU = 26 CPU total. Pinned BLAS keeps it sane.
( export SLURM_CPUS_PER_TASK=6
  Rscript 05g_masl_vs_healthy_dream.R \
    > "${LOGS}/stage_masl_mega.out" 2> "${LOGS}/stage_masl_mega.err"
) &
pid_mega=$!

declare -a pids_loo=()
for cohort in "${COHORTS_MH[@]}"; do
    (
        export CONTRAST="masl_vs_healthy"
        export HELD_OUT="${cohort}"
        export SLURM_CPUS_PER_TASK=4
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

echo "MASL mega rc=${rc_mega}  loo rc=${rc_loo}  $(date)"
if [[ ${rc_mega} -ne 0 || ${rc_loo} -ne 0 ]]; then
    echo "MASL cpu stage FAILED before sweep"
    exit 1
fi

echo "MASL sweep starting $(date)"
( export CONTRAST="masl_vs_healthy" SLURM_CPUS_PER_TASK=2
  Rscript lfc_sweep_loo.R \
    > "${LOGS}/stage_masl_sweep.out" \
    2> "${LOGS}/stage_masl_sweep.err"
)
echo "MASL cpu stage OK $(date)"
