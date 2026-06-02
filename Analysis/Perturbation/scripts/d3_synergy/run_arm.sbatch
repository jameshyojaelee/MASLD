#!/bin/bash
#SBATCH --job-name=D3-runner
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-3        # 4 models
#SBATCH --output=logs/d3_arm_%A_%a.out
#SBATCH --error=logs/d3_arm_%A_%a.err

# D3 Synergy arm — model-runner array.
# Usage: bash ../shared/run_arm_driver.sh d3_synergy
# TIER env var (tier1_pairs / tier2_pairs / tier3_triples / tier4_quadruples / tier5_adaptive)

set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT/Analysis/Perturbation/scripts/d3_synergy"

MODELS=(gears state_pair scgpt_mask tahoe_drug)
MODEL_PRETTY=(GEARS STATE scGPT Tahoe)
MODEL="${MODELS[$SLURM_ARRAY_TASK_ID]}"
PRETTY="${MODEL_PRETTY[$SLURM_ARRAY_TASK_ID]}"
MODALITY="${MODALITY:-zero_shot}"
TIER="${TIER:-tier1_pairs}"
CONTEXT="${CONTEXT:-all}"
EXTRA_FLAGS="${EXTRA_FLAGS:-}"

scontrol update job=$SLURM_JOB_ID name=$PRETTY 2>/dev/null || true

case "$MODEL" in
    gears)      ENV="perturbation_gears" ;;
    state_pair) ENV="perturbation_state" ;;
    scgpt_mask) ENV="perturbation_scgpt" ;;
    tahoe_drug) ENV="perturbation_tahoe" ;;
esac

eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV"
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true
set -u

python "${MODEL}_runner.py" \
    --modality "$MODALITY" \
    --tier "$TIER" \
    --context "$CONTEXT" \
    $EXTRA_FLAGS

echo "[D3/$PRETTY] runner finished"
