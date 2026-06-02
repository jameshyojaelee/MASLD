#!/usr/bin/env bash
# Generic packed LOO runner — one (CONTRAST, MASH_DEF) set with all folds in
# parallel inside a single cpu-partition job. Reads env vars:
#   CONTRAST       mash_vs_masl | mash_vs_healthy | masl_vs_healthy
#   MASH_DEF       borderline_grouped | strict (ignored for masl_vs_healthy)
#   CPU_PER_FOLD   integer (cores given to each parallel R proc)
#   COHORTS_CSV    comma-separated cohort ids to hold out
#   TAG_SUFFIX     optional suffix for log filenames (e.g. "primary" / "strict")
# Pins BLAS to 1 thread per worker so BPPARAM workers actually get their cores.

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/mash_loo_sweep

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# BLAS thread pinning prevents 8-proc × N-thread fanout
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

set -euo pipefail
cd "${SCRIPTS}"

: "${CONTRAST:?CONTRAST env var required}"
: "${CPU_PER_FOLD:?CPU_PER_FOLD env var required}"
: "${COHORTS_CSV:?COHORTS_CSV env var required}"
MASH_DEF=${MASH_DEF:-borderline_grouped}
TAG_SUFFIX=${TAG_SUFFIX:-}

tag="${CONTRAST}"
[[ "${MASH_DEF}" == "strict" ]] && tag="${tag}_strict"

IFS=',' read -r -a cohorts <<< "${COHORTS_CSV}"

echo "LOO set ${tag}  node $(hostname)  cpu=$(nproc)  mem=$(free -g | awk '/^Mem:/ {print $2"G"}')  start $(date)"
echo "  cohorts: ${cohorts[*]}"
echo "  cpu_per_fold: ${CPU_PER_FOLD}  (BLAS pinned to 1)"

declare -a pids=()
for cohort in "${cohorts[@]}"; do
    (
        export CONTRAST="${CONTRAST}"
        export MASH_DEF="${MASH_DEF}"
        export HELD_OUT="${cohort}"
        export SLURM_CPUS_PER_TASK="${CPU_PER_FOLD}"
        Rscript dream_loo_cv_contrasts.R \
            > "${LOGS}/loo_${tag}_${cohort}.out" \
            2> "${LOGS}/loo_${tag}_${cohort}.err"
    ) &
    pids+=($!)
done

set +e
rc_set=0
for p in "${pids[@]}"; do
    wait "${p}" || rc_set=$?
done
set -e

echo "LOO set ${tag} done rc=${rc_set} $(date)"
exit ${rc_set}
