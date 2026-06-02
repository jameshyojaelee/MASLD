#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=env_state
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_state.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Prevent user-site (~/.local/lib/pythonX.Y/site-packages) leak into env Python.
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_state
YAML=envs/perturbation/perturbation_state.yml
export LD_LIBRARY_PATH="$ENV_PATH/lib:${LD_LIBRARY_PATH:-}"

echo "=========================================="
echo "Building $ENV_PATH from $YAML"
echo "Start time: $(date -Iseconds)"
echo "=========================================="

# Remove existing partial env if present
if [ -d "$ENV_PATH" ]; then
  echo "Removing existing partial env at $ENV_PATH"
  rm -rf "$ENV_PATH"
fi

T0=$(date +%s)
# Use solver=libmamba; --yes to skip prompts
micromamba env create -f "$YAML" -p "$ENV_PATH" -y || {
  echo "ERROR: micromamba env create failed for STATE env"
  exit 1
}
T1=$(date +%s)
echo "Solve+install time: $((T1-T0))s"

echo "=== Install arc-state + pertpy via pip with --no-deps ==="
# CRITICAL: --no-deps prevents arc-state from pip-installing torch>=2.7 + nccl-cu13
# which would uninstall conda's torch 2.5.1 and break the env.
# arc-state's non-torch deps (hydra-core, omegaconf, antlr4, cell-eval, cell-load,
# geomloss, pdex, texttable, toml) are not strictly required for import-time smoke tests
# but will need to be added incrementally if downstream code imports them.
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --no-deps arc-state || \
  micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --no-deps git+https://github.com/ArcInstitute/state.git || \
  echo "WARN: arc-state install fallback failed; will need manual follow-up"

# Install arc-state runtime helpers that don't pull torch
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --no-deps \
  hydra-core omegaconf antlr4-python3-runtime cell-eval cell-load geomloss pdex texttable toml || \
  echo "WARN: some arc-state helper deps may be missing"

# pertpy / sc-toolbox / cell2sentence — keep deps (these don't drag a different torch)
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install \
  "pertpy>=0.9" sc-toolbox cell2sentence || \
  echo "WARN: pertpy install partial"

echo "=== Verification ==="
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -c "
import sys, importlib
print('Python:', sys.version)
mods = ['torch','scanpy','anndata','scvi','transformers','lightning','pertpy']
for m in mods:
    try:
        mod = importlib.import_module(m)
        v = getattr(mod, '__version__', '?')
        print(f'  OK {m} {v}')
    except Exception as e:
        print(f'  FAIL {m}: {e}')
try:
    import torch
    print('  CUDA available:', torch.cuda.is_available())
except Exception as e:
    print('  torch import failed:', e)
try:
    import state
    print('  OK state', getattr(state,'__version__','?'))
except Exception as e:
    print(f'  STATE NOT IMPORTABLE (try arc_state): {e}')
    try:
        import arc_state
        print('  OK arc_state', getattr(arc_state,'__version__','?'))
    except Exception as e2:
        print(f'  arc_state also missing: {e2}')
"
echo "End time: $(date -Iseconds)"
echo "DONE: STATE env build"
