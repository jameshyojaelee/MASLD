#!/bin/bash
#SBATCH --job-name=D2-runner
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-3   # 4 models
#SBATCH --output=logs/d2_arm_%A_%a.out
#SBATCH --error=logs/d2_arm_%A_%a.err

# D2 Reversal arm — model-runner array.
# Usage: bash ../shared/run_arm_driver.sh d2_reversal
# REFERENCE env var controls which healthy reference (ref_a/ref_b/ref_c).

set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT/Analysis/Perturbation/scripts/d2_reversal"

MODELS=(scgen cpa state_reversal tahoe_reversal)
MODEL_PRETTY=(scGEN CPA STATE Tahoe)
MODEL="${MODELS[$SLURM_ARRAY_TASK_ID]}"
PRETTY="${MODEL_PRETTY[$SLURM_ARRAY_TASK_ID]}"
MODALITY="${MODALITY:-zero_shot}"
CONTEXT="${CONTEXT:-all}"
REFERENCE="${REFERENCE:-ref_a}"
EXTRA_FLAGS="${EXTRA_FLAGS:-}"

scontrol update job=$SLURM_JOB_ID name=$PRETTY 2>/dev/null || true

case "$MODEL" in
    scgen|cpa)      ENV="perturbation_scgen_cpa" ;;
    state_reversal) ENV="perturbation_state" ;;
    tahoe_reversal) ENV="perturbation_tahoe" ;;
esac

eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV"
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true
set -u

python "${MODEL}_runner.py" \
    --modality "$MODALITY" \
    --context "$CONTEXT" \
    --reference "$REFERENCE" \
    $EXTRA_FLAGS

echo "[D2/$PRETTY] runner finished"
