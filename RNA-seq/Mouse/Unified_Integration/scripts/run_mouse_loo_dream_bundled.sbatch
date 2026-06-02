#!/bin/bash
# run_mouse_loo_dream_bundled.sbatch
# ---------------------------------------------------------------------------
# Bundled version of the per-diet LOO refits: single bigmem allocation runs
# all 5 folds in parallel via background R subprocesses, sharing one
# allocation rather than 5x 210G separate jobs. Each refit needs ~60-80GB,
# so 5 in parallel ~ 300-400GB peak; bigmem allocation comfortably holds it.
#
# Submit:
#   sbatch RNA-seq/Mouse/Unified_Integration/scripts/run_mouse_loo_dream_bundled.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=mouse_loo_bundled
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=40
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/Mouse/Unified_Integration/logs/loo_bundled_%j.log

set -o pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/Mouse/Unified_Integration/logs

# conda activate scripts reference unset env vars (ADDR2LINE etc); avoid set -u here.
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

# Auto-discover diet models from YAML config
YAML_FILE="config/mouse_datasets.yaml"
mapfile -t DIETS < <(python3 -c "
import yaml, sys
with open('$YAML_FILE') as f:
    cfg = yaml.safe_load(f)
diets = cfg.get('diet_models', [])
if not diets:
    print('ERROR: no diet_models in YAML', file=sys.stderr); sys.exit(1)
for d in diets:
    print(d)
")

N_DIETS=${#DIETS[@]}
if [[ $N_DIETS -eq 0 ]]; then
  echo "ERROR: no diet models found in $YAML_FILE" >&2
  exit 1
fi

LOG_DIR="RNA-seq/Mouse/Unified_Integration/logs"

echo "[$(date)] Bundled LOO start ($N_DIETS diets): ${DIETS[*]}"

# Per-fold cores: total CPUs / N folds
CORES_PER_FOLD=$(( ${SLURM_CPUS_PER_TASK:-40} / N_DIETS ))
[[ $CORES_PER_FOLD -lt 1 ]] && CORES_PER_FOLD=1
export OMP_NUM_THREADS=$CORES_PER_FOLD
export OPENBLAS_NUM_THREADS=$CORES_PER_FOLD
export MKL_NUM_THREADS=$CORES_PER_FOLD

declare -a PIDS
for DIET in "${DIETS[@]}"; do
  echo "[$(date)] Launching fold: $DIET"
  LOO_DIET="$DIET" Rscript RNA-seq/Mouse/Unified_Integration/scripts/M_loo_dream_refits.R \
    > "$LOG_DIR/loo_${DIET}.log" 2>&1 &
  PIDS+=($!)
done

# Wait for all folds; collect any failures
RC=0
for i in "${!PIDS[@]}"; do
  if ! wait "${PIDS[$i]}"; then
    echo "[$(date)] FAIL: fold ${DIETS[$i]} (pid ${PIDS[$i]})"
    RC=1
  else
    echo "[$(date)] OK:   fold ${DIETS[$i]}"
  fi
done

echo "[$(date)] Bundled LOO done (rc=$RC)"
exit $RC
