#!/bin/bash
#SBATCH --array=1-30
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=loocv
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/loocv/loocv_extra_%a_%A.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/loocv/loocv_extra_%a_%A.err
# ===========================================================================
# Per-(method x held-out-cohort) LOO-CV array for the 6 NEW DE engines that
# were missing from loocv_perfold.csv:
#   limma_voom, limma_voom_qw, limma_trend, edger_qlf, edger_qlf_robust, edger_lrt
#
# 30 tasks = 6 methods x 5 mega cohorts, ONE method per task (so each writes a
# distinct loocv_<cohort>__<method>.csv via METHOD_TAG in loocv_multimethod.R).
# dream/deseq2/metafor folds already exist on disk (loocv_<cohort>.csv) and are
# NOT recomputed here. The aggregator globs loocv_<cohort>*.csv and rbinds all.
#
# NO %N throttle (repo rule). MAXIMALLY parallel — at 8 cpus/task all 30 fit in
# the idle cpu pool at once. Submit a few of the heaviest folds (edger_*) on
# io/bigmem --qos=interactive too is unnecessary: a single cpu/nslab 30-task
# array starts ~immediately given headroom. The aggregator is dependency-gated.
#
# Modes (dual #!/bin/bash so bash arrays work under sbatch --wrap dash too):
#   bash run_loocv_extra_methods.sh   -> submit the 30-task array + aggregator
#   (SLURM sets SLURM_ARRAY_TASK_ID)  -> run one (method,cohort) fold
# ===========================================================================
set -eo pipefail

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation
LOG_DIR=${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/multimethod_validation/loocv
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
SELF=${SCRIPTS}/run_loocv_extra_methods.sh

# 5 mega cohorts (config/human_datasets.yaml include_in_mega: true).
COHORTS=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)
# 6 new engines (one method per array task; loocv_multimethod tags the output).
METHODS=(limma_voom limma_voom_qw limma_trend edger_qlf edger_qlf_robust edger_lrt)
N_COHORTS=${#COHORTS[@]}   # 5
N_METHODS=${#METHODS[@]}   # 6
N_TASKS=$((N_COHORTS * N_METHODS))   # 30

# ---------------------------------------------------------------------------
# SUBMIT MODE
# ---------------------------------------------------------------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  mkdir -p "${LOG_DIR}"
  echo "=== LOO-CV extra-methods: submitting ${N_TASKS}-task array (6 methods x 5 cohorts) ==="
  echo "Methods: ${METHODS[*]}"
  echo "Cohorts: ${COHORTS[*]}"

  ARRAY_JID=$(env -u SLURM_JOB_ID sbatch --parsable --array=1-${N_TASKS} "$SELF")
  echo "  Array submitted: job ${ARRAY_JID} (tasks 1-${N_TASKS}, cpu/nslab, 8c/64G)"

  AGG_JID=$(env -u SLURM_JOB_ID sbatch --parsable \
      --partition=cpu --qos=nslab \
      --cpus-per-task=4 --mem=16G --time=48:00:00 \
      --job-name=loocv \
      --dependency=afterany:${ARRAY_JID} \
      --output="${LOG_DIR}/loocv_aggregate_extra_%j.out" \
      --error="${LOG_DIR}/loocv_aggregate_extra_%j.err" \
      --wrap "cd ${PROJECT}; export MASLD_PROJECT_ROOT=${PROJECT}; ${RBIN} ${SCRIPTS}/aggregate_loocv_multimethod.R")
  echo "  Aggregation submitted: job ${AGG_JID} (afterany:${ARRAY_JID})"
  echo ""
  echo "Monitor: squeue -u \$(whoami) --name=loocv"
  echo "Map:     task t -> method=METHODS[(t-1)/${N_COHORTS}], cohort=COHORTS[(t-1)%${N_COHORTS}]"
  echo "Results: ${PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/loocv/"
  exit 0
fi

# ---------------------------------------------------------------------------
# ARRAY MODE: one (method, cohort) fold
#   task index t in [1..30]; 0-based i = t-1
#   method index = i / N_COHORTS   (integer div)  -> 6 blocks of 5
#   cohort index = i % N_COHORTS
# ---------------------------------------------------------------------------
cd "${PROJECT}"
set -u
i=$((SLURM_ARRAY_TASK_ID - 1))
M_IDX=$((i / N_COHORTS))
C_IDX=$((i % N_COHORTS))
export MASLD_PROJECT_ROOT="${PROJECT}"
export HELD_OUT="${COHORTS[$C_IDX]}"
export VALIDATION_METHODS="${METHODS[$M_IDX]}"

echo "=== LOO-CV extra fold: METHOD=${VALIDATION_METHODS}  HELD_OUT=${HELD_OUT} ==="
echo "Started: $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}  ARRAY_TASK: ${SLURM_ARRAY_TASK_ID}  (m_idx=${M_IDX}, c_idx=${C_IDX})"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
${RBIN} "${SCRIPTS}/loocv_multimethod.R"
echo "Finished: $(date)"
