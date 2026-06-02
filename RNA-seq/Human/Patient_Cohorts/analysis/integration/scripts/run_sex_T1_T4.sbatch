#!/bin/bash
# run_sex_T1_T4.sbatch
# ---------------------------------------------------------------------------
# Master launcher for Team A sensitivity tests:
#   T1: Power-matched F-control subsampling (B=50)
#   T4: Sample-level bootstrap CI on interaction LFC (B=200; 40 tasks × 5 reps)
#
# This script is meant to be RUN (bash run_sex_T1_T4.sbatch) — it submits the
# four child sbatch jobs (T1 array, T1 aggregate, T4 array, T4 aggregate) with
# afterok dependencies. We submit it via sbatch as a tiny driver job for
# logging consistency.
#
#SBATCH --job-name=sex_T1T4_launch
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/outputs/team_A/A3_launch_slurm_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/outputs/team_A/A3_launch_slurm_%j.err

set -euo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=${PROJ}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOG=${PROJ}/outputs/team_A/A3_launch.log

T1_LOG_DIR=${PROJ}/RNA-seq/results/audit_sensitivity/sex_power_matched/logs
T4_LOG_DIR=${PROJ}/RNA-seq/results/audit_sensitivity/sex_bootstrap/logs
mkdir -p "$T1_LOG_DIR" "$T4_LOG_DIR"

echo "=== A3 launcher: $(date) ===" | tee -a "$LOG"

# ----- T1 array job ---------------------------------------------------------
T1_ARRAY_JID=$(sbatch \
  --parsable \
  --job-name=sex_T1_pmsub \
  --partition=cpu \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task=16 \
  --mem=64G \
  --time=48:00:00 \
  --array=1-50%10 \
  --output=${T1_LOG_DIR}/T1_pm_%A_%a.out \
  --error=${T1_LOG_DIR}/T1_pm_%A_%a.err \
  --wrap "set -eu; \
    eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\"; \
    micromamba activate rnaseq; \
    export R_PARALLEL_SEED=42; \
    export SUBSAMPLE_ITER=\${SLURM_ARRAY_TASK_ID}; \
    export N_PER_CELL_M=52; \
    Rscript ${SCRIPTS}/14.4d_sex_power_matched_T1.R")
echo "T1 array job: ${T1_ARRAY_JID}" | tee -a "$LOG"

# ----- T1 aggregate (depends on T1 array completion) -----------------------
T1_AGG_JID=$(sbatch \
  --parsable \
  --job-name=sex_T1_agg \
  --partition=cpu \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task=2 \
  --mem=16G \
  --time=02:00:00 \
  --dependency=afterok:${T1_ARRAY_JID} \
  --output=${T1_LOG_DIR}/T1_agg_%j.out \
  --error=${T1_LOG_DIR}/T1_agg_%j.err \
  --wrap "set -eu; \
    eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\"; \
    micromamba activate rnaseq; \
    Rscript ${SCRIPTS}/14.4d_aggregate_power_matched_T1.R")
echo "T1 aggregate job: ${T1_AGG_JID} (depends afterok:${T1_ARRAY_JID})" | tee -a "$LOG"

# ----- T4 array job ---------------------------------------------------------
# 40 array tasks × 5 reps each = 200 bootstrap replicates total
T4_ARRAY_JID=$(sbatch \
  --parsable \
  --job-name=sex_T4_boot \
  --partition=cpu \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task=16 \
  --mem=128G \
  --time=72:00:00 \
  --array=1-40%6 \
  --output=${T4_LOG_DIR}/T4_boot_%A_%a.out \
  --error=${T4_LOG_DIR}/T4_boot_%A_%a.err \
  --wrap "set -eu; \
    eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\"; \
    micromamba activate rnaseq; \
    export R_PARALLEL_SEED=42; \
    export BOOT_ITER=\${SLURM_ARRAY_TASK_ID}; \
    export N_BOOT=5; \
    Rscript ${SCRIPTS}/14.4e_sex_bootstrap_T4.R")
echo "T4 array job: ${T4_ARRAY_JID}" | tee -a "$LOG"

# ----- T4 aggregate ---------------------------------------------------------
T4_AGG_JID=$(sbatch \
  --parsable \
  --job-name=sex_T4_agg \
  --partition=cpu \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task=2 \
  --mem=32G \
  --time=02:00:00 \
  --dependency=afterok:${T4_ARRAY_JID} \
  --output=${T4_LOG_DIR}/T4_agg_%j.out \
  --error=${T4_LOG_DIR}/T4_agg_%j.err \
  --wrap "set -eu; \
    eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\"; \
    micromamba activate rnaseq; \
    Rscript ${SCRIPTS}/14.4e_aggregate_bootstrap_T4.R")
echo "T4 aggregate job: ${T4_AGG_JID} (depends afterok:${T4_ARRAY_JID})" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== Submitted jobs: ===" | tee -a "$LOG"
echo "  T1 array     : ${T1_ARRAY_JID}  (sex_power_matched, cpu, 16cpu/64G/48h, array 1-50%10)"  | tee -a "$LOG"
echo "  T1 aggregate : ${T1_AGG_JID}  (afterok:${T1_ARRAY_JID})"                                 | tee -a "$LOG"
echo "  T4 array     : ${T4_ARRAY_JID}  (sex_bootstrap, bigmem, 16cpu/128G/72h, array 1-40%6, 5 reps/task = 200)" | tee -a "$LOG"
echo "  T4 aggregate : ${T4_AGG_JID}  (afterok:${T4_ARRAY_JID})"                                 | tee -a "$LOG"
echo "Done: $(date)" | tee -a "$LOG"
