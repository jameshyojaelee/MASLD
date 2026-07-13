#!/usr/bin/env bash
# ============================================================================
# run_345_346_donorcollapse.sh
#
# Rerun the per-donor LIANA + stage-CCC mixed model with the PSEUDOREPLICATION
# FIX (2026-07-12): pool cells across a true donor's sequencing runs and run
# LIANA ONCE per biological donor, then fit the stage model on true donors.
#
# Writes to NEW locations only -- the run-level (buggy) outputs are preserved:
#   per_donor_lr_donorcollapsed/            (345 per-donor parquets)
#   all_donor_lr_scores_donorcollapsed.tsv.gz  (345b merge)
#   donor_collapsed/stage_lr_lmm_*.tsv      (346 mixed model)
#
# Submit from the login node (only issues sbatch calls):
#   bash Analysis/SingleCell/scripts/run_345_346_donorcollapse.sh
# ============================================================================
set -euo pipefail
BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPTS="$BASE/Analysis/SingleCell/scripts"
ST="$BASE/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
LOGS="$ST/logs"
mkdir -p "$LOGS"
cd "$BASE"

MERGED_TSV="$ST/all_donor_lr_scores_donorcollapsed.tsv.gz"
PER_DONOR_DIR="$ST/per_donor_lr_donorcollapsed"

# 121 true donors / 4 per task = 31 array tasks. --qos=interactive -> %4 ceiling.
# (131 -> 121 on 2026-07-12 when GSE136103's 20 GSMs began pooling to 10 patients.)
N_DONORS=121
BATCH=4
N_TASKS=$(( (N_DONORS + BATCH - 1) / BATCH ))

# ---------- Step 345: per-TRUE-DONOR LIANA (array) -------------------------
JID_345=$(sbatch --parsable \
  --job-name=liana \
  --partition=cpu --qos=interactive \
  --time=48:00:00 --cpus-per-task=8 --mem=64G \
  --array=1-${N_TASKS}%4 \
  --output="$LOGS/345dc_%A_%a.out" --error="$LOGS/345dc_%A_%a.err" \
  --wrap "micromamba run -n rapids_singlecell python $SCRIPTS/345_per_donor_liana.py \
            --array-task-id \$SLURM_ARRAY_TASK_ID --array-batch-size ${BATCH} \
            \$([ \$SLURM_ARRAY_TASK_ID -eq 1 ] && echo --write-counts)")
echo "[submit] 345 donor-collapsed LIANA array (${N_TASKS} tasks) -> $JID_345"

# ---------- Step 345b: merge per-donor parquets ----------------------------
JID_345B=$(sbatch --parsable \
  --dependency=afterok:$JID_345 \
  --job-name=liana \
  --partition=cpu --qos=interactive \
  --time=48:00:00 --cpus-per-task=2 --mem=32G \
  --output="$LOGS/345b_dc_%j.out" --error="$LOGS/345b_dc_%j.err" \
  --wrap "micromamba run -n rapids_singlecell python $SCRIPTS/345b_merge_per_donor_parquets.py \
            --in-dir '$PER_DONOR_DIR' --out '$MERGED_TSV'")
echo "[submit] 345b merge -> $JID_345B (after 345)"

# ---------- Step 346: donor-collapsed mixed-effects (coarse + documented F) -
JID_346=$(sbatch --parsable \
  --dependency=afterok:$JID_345B \
  --job-name=lmer \
  --partition=cpu --qos=interactive \
  --time=48:00:00 --cpus-per-task=16 --mem=128G \
  --output="$LOGS/346_dc_%j.out" --error="$LOGS/346_dc_%j.err" \
  --export=ALL,DONOR_COLLAPSE=TRUE,LR_SCORES_TSV="$MERGED_TSV" \
  --wrap "micromamba run -n rnaseq Rscript $SCRIPTS/346_stage_ccc_mixedmodel.R")
echo "[submit] 346 donor-collapsed mixed model -> $JID_346 (after 345b)"

echo ""
echo "Chain: 345=$JID_345 -> 345b=$JID_345B -> 346=$JID_346"
echo "Monitor: squeue -u \$USER --name=liana,lmer --format='%.10i %.10j %.8T %.10M %R'"
