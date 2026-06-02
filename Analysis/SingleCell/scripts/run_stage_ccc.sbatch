#!/usr/bin/env bash
# ============================================================================
# run_stage_ccc.sbatch
#
# Orchestrate the stage-stratified CCC pipeline (Scripts 343 -> 349).
# Submit with:
#   sbatch run_stage_ccc.sbatch
#
# This script SUBMITS each step as a separate dependent sbatch so the heavy
# LIANA array (345) doesn't block downstream and so each step gets the right
# partition / qos.
#
# Conventions:
#   - cpu --qos=interactive for lmer-heavy R steps (346/347) — bypass 7TB cap
#   - gpu --qos=interactive for the LIANA array (345)
#   - io --qos=interactive for the metadata harvest (343)
#
# Outputs land in:
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/
# ============================================================================
#SBATCH --job-name=stage_ccc_orchestrator
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/orchestrator_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/orchestrator_%j.err

set -euo pipefail
BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPTS="$BASE/Analysis/SingleCell/scripts"
LOGS="$BASE/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs"
mkdir -p "$LOGS"

cd "$BASE"

# ---------- Step 343: harvest documented F-stage --------------------------
JID_343=$(sbatch --parsable \
  --job-name=stage_ccc_343 \
  --partition=cpu --qos=interactive \
  --time=04:00:00 --cpus-per-task=2 --mem=8G \
  --output="$LOGS/343_%j.out" --error="$LOGS/343_%j.err" \
  --wrap "micromamba run -n rapids_singlecell python $SCRIPTS/343_harvest_documented_fstage.py")
echo "[submit] 343 (harvest F-stage) -> $JID_343"

# ---------- Step 345: per-donor LIANA (array, 1 task = 8 donors) ---------
# 269 donors / 8 per task = 34 array tasks (last task processes <8 donors)
# cpu partition, --qos=interactive (4 concurrent ceiling -> array %4).
# Smoke test confirmed CPU 48.5 s/donor -> ~6.5 min per 8-donor task.
JID_345=$(sbatch --parsable \
  --job-name=stage_ccc_345 \
  --partition=cpu --qos=interactive \
  --time=04:00:00 --cpus-per-task=8 --mem=64G \
  --array=1-34%4 \
  --output="$LOGS/345_%A_%a.out" --error="$LOGS/345_%A_%a.err" \
  --wrap "micromamba run -n rapids_singlecell python $SCRIPTS/345_per_donor_liana.py \
            --array-task-id \$SLURM_ARRAY_TASK_ID --array-batch-size 8 \
            \$([ \$SLURM_ARRAY_TASK_ID -eq 1 ] && echo --write-counts)")
echo "[submit] 345 (per-donor LIANA array) -> $JID_345"

# ---------- Step 344: build extended donor metadata (after 343 + 345 counts)
JID_344=$(sbatch --parsable \
  --dependency=afterok:$JID_343:$JID_345 \
  --job-name=stage_ccc_344 \
  --partition=cpu --qos=interactive \
  --time=02:00:00 --cpus-per-task=2 --mem=16G \
  --output="$LOGS/344_%j.out" --error="$LOGS/344_%j.err" \
  --wrap "micromamba run -n rnaseq Rscript $SCRIPTS/344_build_donor_metadata_extended.R")
echo "[submit] 344 (build extended donor meta) -> $JID_344 (after 343 + 345)"

# ---------- Step 346: mixed-effects regression (after 344 + 345) ---------
JID_346=$(sbatch --parsable \
  --dependency=afterok:$JID_344 \
  --job-name=stage_ccc_346 \
  --partition=cpu --qos=interactive \
  --time=24:00:00 --cpus-per-task=16 --mem=128G \
  --output="$LOGS/346_%j.out" --error="$LOGS/346_%j.err" \
  --wrap "micromamba run -n rnaseq Rscript $SCRIPTS/346_stage_ccc_mixedmodel.R")
echo "[submit] 346 (mixed-effects regression) -> $JID_346 (after 344)"

# ---------- Step 347: LOO-dataset replication ----------------------------
JID_347=$(sbatch --parsable \
  --dependency=afterok:$JID_346 \
  --job-name=stage_ccc_347 \
  --partition=cpu --qos=interactive \
  --time=24:00:00 --cpus-per-task=16 --mem=64G \
  --output="$LOGS/347_%j.out" --error="$LOGS/347_%j.err" \
  --wrap "micromamba run -n rnaseq Rscript $SCRIPTS/347_stage_ccc_loo_dataset.R")
echo "[submit] 347 (LOO replication) -> $JID_347 (after 346)"

# ---------- Step 348: bulk concordance ----------------------------------
JID_348=$(sbatch --parsable \
  --dependency=afterok:$JID_346 \
  --job-name=stage_ccc_348 \
  --partition=cpu --qos=interactive \
  --time=04:00:00 --cpus-per-task=4 --mem=32G \
  --output="$LOGS/348_%j.out" --error="$LOGS/348_%j.err" \
  --wrap "micromamba run -n rnaseq Rscript $SCRIPTS/348_stage_ccc_bulk_concordance.R")
echo "[submit] 348 (bulk concordance) -> $JID_348 (after 346)"

# ---------- Step 349: figure -----------------------------------------------
JID_349=$(sbatch --parsable \
  --dependency=afterok:$JID_347:$JID_348 \
  --job-name=stage_ccc_349 \
  --partition=cpu --qos=interactive \
  --time=04:00:00 --cpus-per-task=4 --mem=32G \
  --output="$LOGS/349_%j.out" --error="$LOGS/349_%j.err" \
  --wrap "micromamba run -n rnaseq Rscript $SCRIPTS/349_stage_ccc_figures.R")
echo "[submit] 349 (figures) -> $JID_349 (after 347 + 348)"

echo ""
echo "Submitted chain:"
echo "  343 = $JID_343  (harvest F-stage)"
echo "  345 = $JID_345  (per-donor LIANA array)"
echo "  344 = $JID_344  (extended meta)"
echo "  346 = $JID_346  (mixed-effects)"
echo "  347 = $JID_347  (LOO replication)"
echo "  348 = $JID_348  (bulk concordance)"
echo "  349 = $JID_349  (figures)"
echo ""
echo "Monitor:  squeue -u \$USER --format='%.10i %.30j %.8T %.10M %.6D %R'"
