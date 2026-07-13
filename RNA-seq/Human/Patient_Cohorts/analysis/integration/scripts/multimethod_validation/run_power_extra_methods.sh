#!/usr/bin/env bash
# run_power_extra_methods.sh
# ---------------------------------------------------------------------------
# Add the limma + edgeR engine families to the simulation POWER/FDR grid so the
# multimethod power figure (panelB power_curves / panelB2 fdr_control) includes
# them alongside the already-computed dream / deseq2 / metafor.
#
# Each of the 6 new methods runs as its OWN --array=1-48 on the default nslab
# QOS, with VALIDATION_METHODS=<method> + POWER_OUT_SUFFIX=_<method> so it writes
#   power_grid_<method>_{GRID_TASK}.csv
# on the SAME deterministic seeds (set.seed(1000 + GRID_TASK*100 + rep)) as the
# existing dream/deseq2/metafor shards. Identical sims => own-universe + common-
# universe metrics merge cleanly in aggregate_power_multimethod.R (which globs
# ^power_grid.*_\d+\.csv$, so the new shards are picked up automatically).
#
# nb_params.rds is already cached => no NB-estimation stage. NO %N array throttle
# (repo rule); let QOS / AssocMaxJobsLimit be the ceiling. Methods are spread
# across cpu / bigmem / io partitions so they pack into whichever is open.
# Single-word --job-name=simulation per the task convention.
#
# Usage: bash run_power_extra_methods.sh   (orchestrator only — runs NO compute)
# ---------------------------------------------------------------------------
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation
LOG_DIR=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/power
mkdir -p "${LOG_DIR}"
RPATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

# method -> partition (spread across the 3 open partitions; io mem-capped so keep <=64G)
METHODS=(limma_voom limma_voom_qw limma_trend edger_qlf edger_qlf_robust edger_lrt)
declare -A PART=(
  [limma_voom]=cpu
  [limma_voom_qw]=cpu
  [limma_trend]=bigmem
  [edger_qlf]=bigmem
  [edger_qlf_robust]=io
  [edger_lrt]=cpu
)
# bigmem enforces a per-job MIN memory on nslab too -> give bigmem arms 500G.
declare -A MEM=(
  [limma_voom]=64G
  [limma_voom_qw]=64G
  [limma_trend]=500G
  [edger_qlf]=500G
  [edger_qlf_robust]=64G
  [edger_lrt]=64G
)

echo "=== Power harness: launching 6 extra-method arrays (limma + edgeR) ==="
echo "Log dir: ${LOG_DIR}"
echo ""

ARRAY_JIDS=()
for m in "${METHODS[@]}"; do
  part=${PART[$m]}
  mem=${MEM[$m]}
  JID=$(env -u SLURM_JOB_ID sbatch --parsable \
      --partition=${part} \
      --qos=nslab \
      --cpus-per-task=16 \
      --mem=${mem} \
      --time=90:00:00 \
      --job-name=simulation \
      --array=1-48 \
      --output="${LOG_DIR}/power_${m}_%a_%A.out" \
      --error="${LOG_DIR}/power_${m}_%a_%A.err" \
      --export=ALL \
      --wrap="
set -e
set -u
export MASLD_PROJECT_ROOT=${PROJECT}
export GRID_TASK=\${SLURM_ARRAY_TASK_ID}
export VALIDATION_METHODS=${m}
export POWER_OUT_SUFFIX=_${m}
cd ${PROJECT}
echo \"=== Power sim: method=${m} GRID_TASK=\${GRID_TASK} on \$(hostname) ===\"
echo \"Started: \$(date)\"
echo \"SLURM_JOB_ID: \${SLURM_JOB_ID}  ARRAY_TASK: \${SLURM_ARRAY_TASK_ID}  CPUS: \${SLURM_CPUS_PER_TASK}\"
${RPATH} ${SCRIPTS}/power_multimethod.R
echo \"Finished: \$(date)\"
")
  echo "  [${part} ${mem}] method=${m}  -> array job ${JID} (--array=1-48)"
  ARRAY_JIDS+=("${JID}")
done

# --- Aggregation: afterany on ALL 6 arrays --------------------------------
DEP=$(IFS=:; echo "afterany:${ARRAY_JIDS[*]}")
AGG_JID=$(env -u SLURM_JOB_ID sbatch --parsable \
    --partition=cpu \
    --qos=nslab \
    --cpus-per-task=4 \
    --mem=16G \
    --time=48:00:00 \
    --job-name=simulation \
    --dependency=${DEP} \
    --output="${LOG_DIR}/power_aggregate_extra_%j.out" \
    --error="${LOG_DIR}/power_aggregate_extra_%j.err" \
    --export=ALL \
    --wrap="
set -e
set -u
export MASLD_PROJECT_ROOT=${PROJECT}
cd ${PROJECT}
echo \"=== Power: aggregation (all 9 methods) ===\"
echo \"Started: \$(date)\"
${RPATH} ${SCRIPTS}/aggregate_power_multimethod.R
echo \"Finished: \$(date)\"
")
echo ""
echo "  Aggregation submitted: job ${AGG_JID} (${DEP})"
echo ""
echo "=== All submitted ==="
echo "Method arrays: ${ARRAY_JIDS[*]}"
echo "Aggregator:    ${AGG_JID}"
echo "Monitor: squeue -u \$(whoami) --name=simulation"
echo "Output:  ${PROJECT}/.../multimethod_validation/power/power_summary.csv"
