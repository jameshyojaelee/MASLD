#!/usr/bin/env bash
# Stage 3: 4 LFC sweeps in parallel on cpu partition (lightweight; <1 hour).

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/mash_loo_sweep

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

set -euo pipefail
cd "${SCRIPTS}"

echo "Stage 3 cpu  node $(hostname)  start $(date)"

run_sweep() {
    local contrast="$1" mash_def="$2"
    local tag="${contrast}"
    [[ "${mash_def}" == "strict" ]] && tag="${tag}_strict"
    (
        export CONTRAST="${contrast}"
        export MASH_DEF="${mash_def}"
        export SLURM_CPUS_PER_TASK=2
        Rscript lfc_sweep_loo.R \
            > "${LOGS}/stage3_sweep_${tag}.out" \
            2> "${LOGS}/stage3_sweep_${tag}.err"
    ) &
}

# Canonical MASH definition (2026-05-13 decision): strict NAS >= 5 only.
# Primary (Borderline-grouped) sweep PDFs are deprecated — data CSVs from
# previous runs remain on disk under RNA-seq/results/audit_sensitivity/ for
# sensitivity reference. Uncomment the two primary lines to regenerate.
# run_sweep mash_vs_masl    borderline_grouped
# run_sweep mash_vs_healthy borderline_grouped
run_sweep mash_vs_masl    strict
run_sweep mash_vs_healthy strict

set +e
wait
set -e
echo "Stage 3 cpu OK $(date)"
