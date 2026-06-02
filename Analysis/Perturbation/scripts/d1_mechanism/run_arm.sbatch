#!/bin/bash
#SBATCH --job-name=D1-runner
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-3            # 4 models: state, scgpt, geneformer, tahoe
#SBATCH --output=logs/d1_arm_%A_%a.out
#SBATCH --error=logs/d1_arm_%A_%a.err

# D1 Mechanism arm — model-runner array.
# Usage (preferred — chains consensus automatically):
#   bash ../shared/run_arm_driver.sh d1_mechanism
# Or manually:
#   ARRAY_JOB=$(sbatch --parsable run_arm.sbatch)
#   sbatch --dependency=afterany:$ARRAY_JOB run_arm_consensus.sbatch
#
# P0-X1 fix (2026-05-21): consensus + backtester + atlas_integrator moved to
# run_arm_consensus.sbatch (dependent job) so a single runner failure doesn't
# decapitate the pipeline.
# P1-X2 fix: PYTHONNOUSERSITE=1 + LD_LIBRARY_PATH=$CONDA_PREFIX/lib
# P1-X3 fix: per-task --job-name renamed to model name only (STATE / scGPT / etc.)

# Note: set -u and set -e separately so micromamba activate (which has unbound vars) works
set -o pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
ARM_DIR="$PROJECT_ROOT/Analysis/Perturbation/scripts/d1_mechanism"
cd "$ARM_DIR"

MODELS=(state scgpt geneformer tahoe)
MODEL_PRETTY=(STATE scGPT Geneformer Tahoe)
MODEL="${MODELS[$SLURM_ARRAY_TASK_ID]}"
PRETTY="${MODEL_PRETTY[$SLURM_ARRAY_TASK_ID]}"
MODALITY="${MODALITY:-zero_shot}"
CONTEXT="${CONTEXT:-all}"
EXTRA_FLAGS="${EXTRA_FLAGS:-}"

# Rename job to model-name only (per user pref feedback_slurm_job_naming.md)
scontrol update job=$SLURM_JOB_ID name=$PRETTY 2>/dev/null || true

case "$MODEL" in
    state)      ENV="perturbation_state" ;;
    scgpt)      ENV="perturbation_scgpt" ;;
    geneformer) ENV="perturbation_state" ;;   # Geneformer shares state env
    tahoe)      ENV="perturbation_tahoe" ;;
    *) echo "Unknown model $MODEL"; exit 1 ;;
esac

echo "[D1/$PRETTY] activating $ENV"
eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV"

# P1-X2: env hygiene
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
unset PYTHONPATH || true

set -u  # safe to enable now that env activate is done

python "${MODEL}_runner.py" \
    --modality "$MODALITY" \
    --context "$CONTEXT" \
    $EXTRA_FLAGS

echo "[D1/$PRETTY] runner finished"
