#!/bin/bash
#SBATCH --job-name=smoke-st-v1
#SBATCH --partition=io
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/state_train_v1/logs/smoke_st_v1_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/state_train_v1/logs/smoke_st_v1_%j.err

# CPU smoke test of `state tx train` configuration — runs 2 training steps
# end-to-end so we catch any Hydra / cell-load / config bugs BEFORE the multi-day
# real training job (potentially queued for hours) starts and fails.
set -eo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
RUN_ROOT="$PROJECT_ROOT/data/perturbation/state_train_v1"

eval "$(micromamba shell hook --shell bash)"
micromamba activate "$PROJECT_ROOT/.envs/perturbation_state"

export PYTHONNOUSERSITE=1
export PYTHONPATH="$PROJECT_ROOT/Analysis/Perturbation/scripts/shared/env_setup/vendor"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export WANDB_MODE=disabled
set -u

echo "[smoke] start: $(date -Is) on $(hostname)"

# Use a throwaway run name so the real run dir isn't polluted
python "$PROJECT_ROOT/Analysis/Perturbation/scripts/shared/env_setup/train_st_head_v1.py" \
    --toml "$RUN_ROOT/toml/replogle_norman_v1.toml" \
    --name "smoke_test_$(date +%s)" \
    --output-dir "$RUN_ROOT/checkpoints/smoke" \
    --embed-key X_hvg \
    --output-space gene \
    --num-workers 2 \
    --batch-size 16 \
    --max-steps 2 \
    --val-freq 2 \
    --ckpt-every-n-steps 2 \
    --overwrite

echo "[smoke] complete: $(date -Is)"
