#!/bin/bash
#SBATCH --job-name=train-st-v1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --time=72:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=12
# NB: l40s (sm_89, 48 GB VRAM) is the right balance for this env's torch
# 2.5.1+cu118 build — b6k (Blackwell) needs newer torch and likely trips
# task-#21 NCCL/CXXABI issues. Mem=96G is comfortable for two ~1-2 GB h5ad
# corpora + 12 dataloader workers; the prior 256 G ask was rejected as
# unsatisfiable when the node mix is saturated. gpu_interactive partition
# is srun-only — can't sbatch there.
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/state_train_v1/logs/train_st_v1_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/state_train_v1/logs/train_st_v1_%j.err

# Train a new STATE Transition (ST) head with expanded Perturb-seq vocabulary,
# Combines:
#   * Replogle 2022 K562 essential (1,093 single-gene KD perturbations, ~163k cells)
#   * Norman 2019 A549/K562 (282 perts, mix of singles + combinatorial, ~91k cells)
#
# Target vocab union ≈ 1,300 perturbations (vs ~2,024 in the prebuilt HepG2 ST
# head). The marginal gain is the combinatorial Norman tokens, which the
# prebuilt heads don't carry. v2 will add Adamson + full genome-wide Replogle.
#
# Architecture: matches the published st-se-replogle-full HepG2 0.99 config
# (cell_set_len=64, hidden_dim=328, energy loss). Trains directly from HVG
# expression (no SE-600M re-embedding upfront — saves 1-2 GPU-hours of prep).
#
# Submit:
#   # 1. Build the corpus (CPU, ~15 min)
#   CORPUS_JOB=$(sbatch --parsable Analysis/Perturbation/scripts/shared/env_setup/assemble_st_train_corpus_v1.sbatch)
#   # 2. Chain training (GPU, multi-day)
#   sbatch --dependency=afterok:$CORPUS_JOB Analysis/Perturbation/scripts/shared/env_setup/train_st_head_v1.sbatch
set -eo pipefail   # `set -u` deferred — micromamba activate scripts reference
                   # unbound vars (MKL_INTERFACE_LAYER) and trip nounset early.

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
RUN_ROOT="$PROJECT_ROOT/data/perturbation/state_train_v1"
mkdir -p "$RUN_ROOT/logs"

eval "$(micromamba shell hook --shell bash)"
micromamba activate "$PROJECT_ROOT/.envs/perturbation_state"

# Strict env-isolation per project policy:
#   * PYTHONNOUSERSITE blocks ~/.local pollution (e.g. torchvision ABI mismatch)
#   * PYTHONPATH points to the vendored cell-load (Arc package not on PyPI in
#     a form compatible with this env; vendoring keeps us from `pip install`
#     into perturbation_state per the project moratorium)
#   * LD_LIBRARY_PATH bias-toward-env for CUDA/NCCL libs
export PYTHONNOUSERSITE=1
export PYTHONPATH="$PROJECT_ROOT/Analysis/Perturbation/scripts/shared/env_setup/vendor"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
set -u   # safe to enable now that env is active

# Disable W&B (no creds wired; project policy: log to local stdout + lightning
# CSV/JSON loggers, which `state tx train` falls back to with use_wandb=false)
export WANDB_MODE=disabled

echo "[train-st-v1] start: $(date -Is)"
echo "[train-st-v1] node: $(hostname); cuda devices: ${CUDA_VISIBLE_DEVICES:-unset}"
echo "[train-st-v1] gpu info:"
nvidia-smi --query-gpu=name,driver_version,memory.total,memory.free --format=csv,noheader || true
python -c "import torch; print('[train-st-v1] torch=%s cuda=%s avail=%s ndevice=%d' % (torch.__version__, torch.version.cuda, torch.cuda.is_available(), torch.cuda.device_count()))"

# Verify the corpus the prep job emitted
ls -1 "$RUN_ROOT/datasets/replogle/" "$RUN_ROOT/datasets/norman/" || {
    echo "[train-st-v1] FATAL: training corpus not assembled — run assemble_st_train_corpus_v1.sbatch first." >&2
    exit 2
}

# Run name + output dir
RUN_NAME="${RUN_NAME:-masld_replogle_norman_v1}"
OUTPUT_DIR="$RUN_ROOT/checkpoints"
mkdir -p "$OUTPUT_DIR"

# Hyperparameters — drawn from the published st-se-replogle-full HepG2 0.99
# config, with these overrides:
#   * embed_key=X_hvg (no SE-600M re-embedding upfront — train directly on
#     HVG expression; SE-based variant available as v2)
#   * cell_type_key=cell_type and pert_col=gene (cell-load defaults)
#   * control_pert=non-targeting
#   * batch_col=gem_group
#   * max_steps=80000 (matches HepG2 0.99; ~36h on L40S at ~1.5 steps/sec)
python "$PROJECT_ROOT/Analysis/Perturbation/scripts/shared/env_setup/train_st_head_v1.py" \
    --toml "$RUN_ROOT/toml/replogle_norman_v1.toml" \
    --name "$RUN_NAME" \
    --output-dir "$OUTPUT_DIR" \
    --embed-key X_hvg \
    --output-space gene \
    --num-workers 8 \
    --batch-size 64 \
    --max-steps 80000 \
    --val-freq 2000 \
    --ckpt-every-n-steps 4000

echo "[train-st-v1] complete: $(date -Is)"
ls -lh "$OUTPUT_DIR/$RUN_NAME/checkpoints/" || true
