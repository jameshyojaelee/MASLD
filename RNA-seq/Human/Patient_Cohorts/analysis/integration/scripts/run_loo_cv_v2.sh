#!/usr/bin/env bash
# run_loo_cv_v2.sh
# ---------------------------------------------------------------------------
# Launch LOO-CV v2: 5 independent folds (one per mega cohort) + aggregation.
# Each fold recomputes filterByExpr + calcNormFactors (Fix A) and tests
# held-out replication (Fix B).
#
# Usage: bash run_loo_cv_v2.sh
# ---------------------------------------------------------------------------
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOG_DIR=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/loo_cv_v2
mkdir -p "${LOG_DIR}"

# 5 mega cohorts (from config/human_datasets.yaml include_in_mega: true)
COHORTS=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)

echo "=== LOO-CV v2: submitting ${#COHORTS[@]} folds ==="
echo "Cohorts: ${COHORTS[*]}"
echo "Log dir: ${LOG_DIR}"
echo ""

FOLD_JIDS=()

for cohort in "${COHORTS[@]}"; do
    JID=$(sbatch \
        --partition=cpu \
        --cpus-per-task=16 \
        --mem=64G \
        --time=48:00:00 \
        --qos=nslab \
        --job-name=dream-loocv \
        --output="${LOG_DIR}/loo_v2_${cohort}_%j.out" \
        --error="${LOG_DIR}/loo_v2_${cohort}_%j.err" \
        --export=ALL,HELD_OUT="${cohort}" \
        --wrap="
set -eo pipefail
eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\"
micromamba activate rnaseq
set -u
cd ${SCRIPTS}
echo \"=== LOO-CV v2 fold: HELD_OUT=${cohort} ===\"
echo \"Started: \$(date)\"
echo \"SLURM_JOB_ID: \${SLURM_JOB_ID}\"
echo \"SLURM_CPUS_PER_TASK: \${SLURM_CPUS_PER_TASK}\"
Rscript dream_loo_cv_v2.R
echo \"Finished: \$(date)\"
" | awk '{print $4}')

    echo "  ${cohort}: submitted job ${JID}"
    FOLD_JIDS+=("${JID}")
done

# Chain aggregation after all folds complete
DEP_STR=$(IFS=:; echo "${FOLD_JIDS[*]}")
AGG_JID=$(sbatch \
    --partition=cpu \
    --cpus-per-task=4 \
    --mem=16G \
    --time=1:00:00 \
    --qos=nslab \
    --job-name=dream-loocv-agg \
    --dependency=afterok:${DEP_STR} \
    --output="${LOG_DIR}/loo_v2_aggregate_%j.out" \
    --error="${LOG_DIR}/loo_v2_aggregate_%j.err" \
    --wrap="
set -eo pipefail
eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\"
micromamba activate rnaseq
set -u
cd ${SCRIPTS}
echo \"=== LOO-CV v2 aggregation ===\"
echo \"Started: \$(date)\"
Rscript aggregate_loo_cv_v2.R
echo \"Finished: \$(date)\"
" | awk '{print $4}')

echo ""
echo "  Aggregation: submitted job ${AGG_JID} (depends on ${DEP_STR})"
echo ""
echo "=== All jobs submitted ==="
echo "Monitor: squeue -u \$(whoami) --name=dream-loocv,dream-loocv-agg"
echo "Results: ${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv_v2/"
