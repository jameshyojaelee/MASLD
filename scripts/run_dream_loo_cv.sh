#!/bin/bash
# run_dream_loo_cv.sh
# Submit 8 parallel SLURM jobs for dream leave-one-study-out cross-validation.
# Each job holds out one cohort and re-runs dream on the remaining cohorts.
# Usage: bash scripts/run_dream_loo_cv.sh

set -euo pipefail

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
LOG_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs"
mkdir -p "$LOG_DIR"

# 8 cohorts used in dream mega-analysis (excludes GSE167523 and PRJNA512027)
COHORTS=(
  GSE126848
  GSE130970
  GSE135251
  GSE162694
  GSE174478
  GSE193066
  GSE213621
  GSE240729
)

echo "Submitting 8 LOO-CV dream jobs..."
JOB_IDS=""

WRAPPER="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/run_single_loo.sh"

for study in "${COHORTS[@]}"; do
  JOB_ID=$(sbatch --parsable \
    --job-name="loo_${study}" \
    --output="${LOG_DIR}/loo_cv_${study}_%j.out" \
    --error="${LOG_DIR}/loo_cv_${study}_%j.err" \
    --export=ALL,HELD_OUT="${study}" \
    "${WRAPPER}")

  echo "  ${study}: Job ${JOB_ID}"
  JOB_IDS="${JOB_IDS}${JOB_ID}:"
done

# Remove trailing colon
JOB_IDS="${JOB_IDS%:}"
echo ""
echo "All 8 LOO-CV jobs submitted."
echo "Job IDs: ${JOB_IDS}"
echo ""

# Auto-chain post-processing after all LOO jobs complete
POST_JOB=$(sbatch --parsable \
  --dependency="afterok:${JOB_IDS}" \
  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/run_loo_cv_post.sh)
echo "Post-processing job ${POST_JOB} queued (depends on all 8 LOO jobs)"
echo "Monitor with: squeue -u \$USER --name='loo_*'"
