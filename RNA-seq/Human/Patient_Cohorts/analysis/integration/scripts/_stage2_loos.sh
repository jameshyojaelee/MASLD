#!/usr/bin/env bash
# Stage 2: 26 LOO folds (4 contrast x mash_def sets) packed inside one
# bigmem node. Within each set, all folds run in parallel.

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/mash_loo_sweep

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -euo pipefail
cd "${SCRIPTS}"

echo "Stage 2 node $(hostname)  cpu=$(nproc)  mem=$(free -g | awk '/^Mem:/ {print $2"G"}')  start $(date)"

COHORTS_MM=(GSE126848 GSE130970 GSE135251 GSE162694 GSE167523 GSE174478 GSE193066)
COHORTS_MH=(GSE126848 GSE130970 GSE135251 GSE162694)

run_loo_set() {
    local contrast="$1" mash_def="$2" cpu_per_fold="$3"
    shift 3
    local tag="${contrast}"
    [[ "${mash_def}" == "strict" ]] && tag="${tag}_strict"
    echo ""
    echo "==== LOO ${tag} (cpu_per_fold=${cpu_per_fold}, folds=$#) start $(date) ===="
    local pids=()
    local rc_set=0
    for cohort in "$@"; do
        (
            export CONTRAST="${contrast}"
            export MASH_DEF="${mash_def}"
            export HELD_OUT="${cohort}"
            export SLURM_CPUS_PER_TASK="${cpu_per_fold}"
            Rscript dream_loo_cv_contrasts.R \
                > "${LOGS}/stage2_${tag}_${cohort}.out" \
                2> "${LOGS}/stage2_${tag}_${cohort}.err"
        ) &
        pids+=($!)
    done
    set +e
    for p in "${pids[@]}"; do
        wait "${p}" || rc_set=$?
    done
    set -e
    echo "==== LOO ${tag} done rc=${rc_set} $(date) ===="
    return ${rc_set}
}

# 8-fold sets: 8 procs * 8 CPU = 64 CPU concurrent
run_loo_set mash_vs_masl    borderline_grouped 8  "${COHORTS_MM[@]}"
run_loo_set mash_vs_masl    strict             8  "${COHORTS_MM[@]}"

# 5-fold sets: 5 procs * 14 CPU = 70 CPU concurrent
run_loo_set mash_vs_healthy borderline_grouped 14 "${COHORTS_MH[@]}"
run_loo_set mash_vs_healthy strict             14 "${COHORTS_MH[@]}"

echo "Stage 2 OK $(date)"
