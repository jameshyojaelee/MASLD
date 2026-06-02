#!/bin/bash
#SBATCH --job-name=fig_regen_lfc05
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=24
#SBATCH --mem=300G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_regen_lfc05_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig_regen_lfc05_%j.err

# Re-render Fig 1 + figS02 + figS_threshold_landscape post LFC=0.5 migration.
# Two phases of within-node parallelism, sharing one bigmem allocation:
#   Phase A: producer scripts (130, 115) emit threshold-tagged contrast tables.
#   Phase B: consumer figure scripts (fig1, figS02, figS_threshold_landscape)
#            each render in parallel since they have no inter-dependencies.
# fig1_atlas_overview_v2 is technically independent of A but kept inside Phase B
# for simplicity (a few extra seconds).

set -o pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
mkdir -p scripts/figures/logs

# conda activate references unset env vars; defer set -u until after.
source "$HOME/.bashrc"
micromamba activate rnaseq
set -u

# Per-subprocess core caps so 5 parallel R sessions don't oversubscribe the 24-core allocation.
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

LOG=scripts/figures/logs

run_step () {
  local label="$1" script="$2"
  echo "[$(date)] Launching: $label ($script)"
  Rscript "$script" > "$LOG/regen_${label}.log" 2>&1
  echo "[$(date)] Done: $label rc=$?"
}

echo "=== Phase A: producer scripts (parallel) ==="
run_step "130_progression" "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/130_progression_contrasts_dream.R" &
PID_130=$!
run_step "115_transition" "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/115_transition_programs.R" &
PID_115=$!
wait $PID_130 || echo "[WARN] 130 failed"
wait $PID_115 || echo "[WARN] 115 failed"

echo "=== Phase B: consumer figure scripts (parallel) ==="
run_step "fig1" "scripts/figures/fig1_atlas_overview_v2.R" &
PID_F1=$!
run_step "figS02" "scripts/figures/figS02_progression_contrasts.R" &
PID_S02=$!
run_step "figS_landscape" "scripts/figures/figS_deg_threshold_landscape.R" &
PID_LAND=$!
wait $PID_F1   || echo "[WARN] fig1 failed"
wait $PID_S02  || echo "[WARN] figS02 failed"
wait $PID_LAND || echo "[WARN] figS_landscape failed"

echo "[$(date)] All fig_regen tasks done."
