#!/bin/bash
#SBATCH --job-name=D5-runner
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-2       # 3 models
#SBATCH --output=logs/d5_arm_%A_%a.out
#SBATCH --error=logs/d5_arm_%A_%a.err

# D5 Mouse cross-species arm — model-runner array.
# Usage: bash ../shared/run_arm_driver.sh d5_mouse

set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$PROJECT_ROOT/Analysis/Perturbation/scripts/d5_mouse"

MODELS=(geneformer_mouse state_mouse scgen_mouse)
MODEL_PRETTY=(Geneformer STATE scGEN)
MODEL="${MODELS[$SLURM_ARRAY_TASK_ID]}"
PRETTY="${MODEL_PRETTY[$SLURM_ARRAY_TASK_ID]}"
MODALITY="${MODALITY:-zero_shot}"
CONTEXT="${CONTEXT:-all}"
EXTRA_FLAGS="${EXTRA_FLAGS:-}"

scontrol update job=$SLURM_JOB_ID name=$PRETTY 2>/dev/null || true

case "$MODEL" in
    geneformer_mouse|state_mouse) ENV="perturbation_state" ;;
    scgen_mouse)                  ENV="perturbation_scgen_cpa" ;;
esac

eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV"
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true
set -u

python "${MODEL}_runner.py" --modality "$MODALITY" --context "$CONTEXT" $EXTRA_FLAGS

echo "[D5/$PRETTY] runner finished"
