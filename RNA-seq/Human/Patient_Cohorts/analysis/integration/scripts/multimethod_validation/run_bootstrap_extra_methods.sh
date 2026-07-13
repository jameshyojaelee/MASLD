#!/usr/bin/env bash
# run_bootstrap_extra_methods.sh
# ---------------------------------------------------------------------------
# Extend the Pillar-A bootstrap (CPSS selection-stability) harness to the limma
# + edgeR engine families so the stability figure (bootstrap_stability_counts /
# selection_frequency / method_concordance) includes them alongside the already-
# computed dream / deseq2 / metafor.
#
# The original run only computed VALIDATION_METHODS=dream,deseq2 (iter_{N}.csv)
# plus a separate metafor arm (iter_metafor_{N}.csv). Each new method runs as its
# OWN --array=1-500 with VALIDATION_METHODS=<method> + BOOTSTRAP_OUT_SUFFIX=_<method>
# so it writes
#   iter_<method>_{ITER}.csv
# on the SAME deterministic per-iteration subsample seed (set.seed(42+ITER*7919)
# in bootstrap_multimethod.R) as the existing dream/deseq2/metafor iters. Identical
# splits => the per-method selection frequencies are directly comparable, and
# aggregate_bootstrap_multimethod.R (globs ^iter.*_\d+\.csv$) picks the new shards
# up automatically and merges all 9 methods.
#
# Resumable: bootstrap_multimethod.R skips any iter file that already exists.
# NO %N array throttle (repo rule); let QOS / AssocMaxJobsLimit be the ceiling.
# Single-word --job-name=multimethod (matches the existing bootstrap study).
# nslab QOS (NOT interactive — that has a 4-job cap). Spread cpu/io (io fits 64G);
# NEVER bigmem (its 500G floor balloons the nslab 7TB memory cap).
#
# Usage: bash run_bootstrap_extra_methods.sh   (orchestrator only — runs NO compute)
# ---------------------------------------------------------------------------
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation
LOG_DIR=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs
mkdir -p "${LOG_DIR}"
RPATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

METHODS=(limma_voom limma_voom_qw limma_trend edger_qlf edger_qlf_robust edger_lrt)
declare -A PART=(
  [limma_voom]=cpu
  [limma_voom_qw]=cpu
  [limma_trend]=cpu
  [edger_qlf]=io
  [edger_qlf_robust]=io
  [edger_lrt]=io
)

echo "=== Bootstrap CPSS: launching 6 extra-method arrays (limma + edgeR) ==="
echo "Log dir: ${LOG_DIR}"
echo ""

ARRAY_JIDS=()
for m in "${METHODS[@]}"; do
  part=${PART[$m]}
  JID=$(env -u SLURM_JOB_ID -u SLURM_ARRAY_TASK_ID sbatch --parsable \
      --partition=${part} \
      --qos=nslab \
      --cpus-per-task=16 \
      --mem=64G \
      --time=48:00:00 \
      --job-name=multimethod \
      --array=1-500 \
      --output="${LOG_DIR}/bootstrap_${m}_%a_%A.out" \
      --error="${LOG_DIR}/bootstrap_${m}_%a_%A.err" \
      --export=ALL \
      --wrap="
set -e
set -u
export MASLD_PROJECT_ROOT=${PROJECT}
export ITER=\${SLURM_ARRAY_TASK_ID}
export VALIDATION_METHODS=${m}
export BOOTSTRAP_OUT_SUFFIX=_${m}
cd ${PROJECT}
${RPATH} ${SCRIPTS}/bootstrap_multimethod.R
")
  echo "  [${part}] method=${m}  -> array job ${JID} (--array=1-500, suffix=_${m})"
  ARRAY_JIDS+=("${JID}")
done

# --- Aggregation: afterany on ALL 6 arrays --------------------------------
DEP=$(IFS=:; echo "afterany:${ARRAY_JIDS[*]}")
AGG_JID=$(env -u SLURM_JOB_ID -u SLURM_ARRAY_TASK_ID sbatch --parsable \
    --partition=cpu --qos=nslab \
    --cpus-per-task=4 --mem=32G --time=48:00:00 \
    --job-name=multimethod \
    --dependency=${DEP} \
    --output="${LOG_DIR}/bootstrap_aggregate_extra_%j.out" \
    --error="${LOG_DIR}/bootstrap_aggregate_extra_%j.err" \
    --export=ALL \
    --wrap="
set -e
set -u
export MASLD_PROJECT_ROOT=${PROJECT}
cd ${PROJECT}
${RPATH} ${SCRIPTS}/aggregate_bootstrap_multimethod.R
")
echo ""
echo "  Aggregation submitted: job ${AGG_JID} (${DEP})"
echo "=== All submitted ==="
echo "Method arrays: ${ARRAY_JIDS[*]}"
echo "Aggregator:    ${AGG_JID}"
