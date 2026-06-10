#!/usr/bin/env bash
# run_lvqw_loo_cv_C2.sh — 5 LVQW LOO-CV (C2) folds + aggregation.
# One fold per mega cohort; each refits the LVQW engine (C2 design) on the
# other 4 cohorts and recovers vs the full C2 canonical.
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOG_DIR=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/loo_cv_C2
mkdir -p "${LOG_DIR}"

COHORTS=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)
PART=${PART:-cpu}

echo "=== LVQW LOO-CV (C2): submitting ${#COHORTS[@]} folds on partition ${PART} ==="
FOLD_JIDS=()
for cohort in "${COHORTS[@]}"; do
    JID=$(env -u SLURM_JOB_ID sbatch \
        --partition=${PART} \
        --qos=interactive \
        --cpus-per-task=16 \
        --mem=64G \
        --time=48:00:00 \
        --job-name=loocv \
        --output="${LOG_DIR}/loo_C2_${cohort}_%j.out" \
        --error="${LOG_DIR}/loo_C2_${cohort}_%j.err" \
        --export=ALL,HELD_OUT="${cohort}" \
        --wrap="bash -c '
eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\"
micromamba activate rnaseq
cd ${SCRIPTS}
echo \"=== LVQW LOO-CV (C2) fold: HELD_OUT=${cohort} ===\"
echo \"Started: \$(date)  JOB=\${SLURM_JOB_ID}  CPUS=\${SLURM_CPUS_PER_TASK}\"
Rscript lvqw_loo_cv_C2.R
'" | awk '{print $4}')
    echo "  ${cohort}: job ${JID}"
    FOLD_JIDS+=("${JID}")
done

DEP_STR=$(IFS=:; echo "${FOLD_JIDS[*]}")
AGG_JID=$(env -u SLURM_JOB_ID sbatch \
    --partition=${PART} \
    --qos=interactive \
    --cpus-per-task=4 \
    --mem=16G \
    --time=4:00:00 \
    --job-name=loocv \
    --dependency=afterok:${DEP_STR} \
    --output="${LOG_DIR}/loo_C2_aggregate_%j.out" \
    --error="${LOG_DIR}/loo_C2_aggregate_%j.err" \
    --wrap="bash -c '
eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\"
micromamba activate rnaseq
cd ${SCRIPTS}
echo \"=== LVQW LOO-CV (C2) aggregation ===  \$(date)\"
Rscript aggregate_lvqw_loo_cv_C2.R
'" | awk '{print $4}')

echo ""
echo "  Aggregation: job ${AGG_JID} (afterok:${DEP_STR})"
echo "Monitor: squeue -u \$(whoami) --name=loocv"
echo "Results: ${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv_C2/"
echo "${AGG_JID}" > "${LOG_DIR}/.agg_jid"
