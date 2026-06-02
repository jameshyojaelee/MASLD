#!/bin/bash
#SBATCH --job-name=D4-runner
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-2     # 3 propagation methods
#SBATCH --output=logs/d4_arm_%A_%a.out
#SBATCH --error=logs/d4_arm_%A_%a.err

# D4 Circuit arm — propagation-method array.
# Usage: bash ../shared/run_arm_driver.sh d4_circuit
# Note: lr_curator now runs in run_arm_consensus.sbatch (P1-D4-5 fix: was racing
# with niches/commot tasks under old array-task-0 special-case).

set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT/Analysis/Perturbation/scripts/d4_circuit"

METHODS=(state_2stage niches commot)
METHOD_PRETTY=(STATE NICHES COMMOT)
METHOD="${METHODS[$SLURM_ARRAY_TASK_ID]}"
PRETTY="${METHOD_PRETTY[$SLURM_ARRAY_TASK_ID]}"
MODALITY="${MODALITY:-zero_shot}"
CONTEXT="${CONTEXT:-all}"
EXTRA_FLAGS="${EXTRA_FLAGS:-}"

scontrol update job=$SLURM_JOB_ID name=$PRETTY 2>/dev/null || true

case "$METHOD" in
    state_2stage) ENV="perturbation_state" ;;
    niches)       ENV="perturbation_ccc" ;;
    commot)       ENV="perturbation_ccc" ;;
esac

eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV"
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true
set -u

if [ "$METHOD" = "niches" ]; then
    # Sanity: NICHES R install (task #20) must be done
    if ! Rscript -e 'library(NICHES)' 2>/dev/null; then
        echo "[D4/NICHES] ERROR: NICHES not installed in perturbation_ccc env. Run task #20 first."
        exit 1
    fi
    Rscript niches_runner.R --modality "$MODALITY" --context "$CONTEXT" $EXTRA_FLAGS
else
    python "${METHOD}_runner.py" --modality "$MODALITY" --context "$CONTEXT" $EXTRA_FLAGS
fi

echo "[D4/$PRETTY] runner finished"
