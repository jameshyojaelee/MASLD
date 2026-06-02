#!/bin/bash
# launch_phase1_5cohort.sh — submit Phase 1 chain (5-cohort dream rebuild).
#
# Chain:
#   A:  05_dream_mega_analysis.R       (bigmem, ~30-60 min)
#   B:  05b_ashr_shrinkage.R           (after A, cpu)
#   C:  dream_loo_cv.R x 5 cohorts     (after A, bigmem, parallel)
#   D:  aggregate_loo_cv.R             (after all C)
#   E:  fig1g_loo_lfc_stability.R      (after D)
#
# All sbatch headers are in scripts/sbatch/_phase1_*.sbatch.
# Cohort set is read from config/human_datasets.yaml at submit time.

set -eo pipefail

REPO_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_DIR=$REPO_ROOT/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
SBATCH_DIR=$SCRIPT_DIR/sbatch

cd "$SCRIPT_DIR"
mkdir -p logs

# 5 mega cohorts (yaml include_in_mega=true). Hardcoded here so we can submit
# 5 LOO jobs from bash; the underlying R scripts re-read the yaml.
COHORTS=(GSE126848 GSE130970 GSE135251 GSE162694 GSE213621)

echo "=== Phase 1 chain submit ==="
echo "Repo:    $REPO_ROOT"
echo "Cohorts: ${COHORTS[*]}"
echo

# --- A: dream mega-analysis ----------------------------------------------
A=$(sbatch --parsable "$SBATCH_DIR/_phase1_05_dream.sbatch")
echo "A  dream:           jobid=$A"

# --- B: ashr (after A) ---------------------------------------------------
B=$(sbatch --parsable --dependency=afterok:$A "$SBATCH_DIR/_phase1_05b_ashr.sbatch")
echo "B  ashr:            jobid=$B  (after $A)"

# --- C1..C5: LOO (after A, parallel) ------------------------------------
C_IDS=()
for COH in "${COHORTS[@]}"; do
  CID=$(sbatch --parsable \
        --job-name=loo_$COH \
        --dependency=afterok:$A \
        --export=ALL,HELD_OUT=$COH \
        "$SBATCH_DIR/_phase1_loo_cv.sbatch")
  echo "C  loo $COH:    jobid=$CID  (after $A)"
  C_IDS+=("$CID")
done
C_DEP=$(IFS=:; echo "${C_IDS[*]}")

# --- D: aggregate_loo_cv (after all C) ----------------------------------
D=$(sbatch --parsable --dependency=afterok:$C_DEP "$SBATCH_DIR/_phase1_aggregate_loo.sbatch")
echo "D  agg_loo:         jobid=$D  (after ${C_IDS[*]})"

# --- E: fig1g (after D) -------------------------------------------------
E=$(sbatch --parsable --dependency=afterok:$D "$SBATCH_DIR/_phase1_fig1g.sbatch")
echo "E  fig1g:           jobid=$E  (after $D)"

echo
echo "=== Submitted job IDs ==="
echo "A=$A"
echo "B=$B"
for i in "${!C_IDS[@]}"; do
  echo "C${i}_${COHORTS[$i]}=${C_IDS[$i]}"
done
echo "D=$D"
echo "E=$E"
echo
echo "All ids (squeue): $A,$B,$(IFS=,; echo "${C_IDS[*]}"),$D,$E"
