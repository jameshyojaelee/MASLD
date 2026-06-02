#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=env_scgpt
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_scgpt.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Prevent user-site (~/.local/lib/pythonX.Y/site-packages) leak into env Python.
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_scgpt
YAML=envs/perturbation/perturbation_scgpt.yml
export LD_LIBRARY_PATH="$ENV_PATH/lib:${LD_LIBRARY_PATH:-}"

echo "=========================================="
echo "Building $ENV_PATH from $YAML"
echo "Start time: $(date -Iseconds)"
echo "NOTE: scGPT may need torch 1.13 + CUDA 11.7 + flash-attn<1.0.5"
echo "      First attempt: as-written torch>=2.0 + flash-attn>=2.5"
echo "      If failure: fallback YAML triggered"
echo "=========================================="

if [ -d "$ENV_PATH" ]; then
  echo "Removing existing partial env at $ENV_PATH"
  rm -rf "$ENV_PATH"
fi

T0=$(date +%s)
micromamba env create -f "$YAML" -p "$ENV_PATH" -y || {
  echo "ERROR: micromamba env create failed for scGPT env"
  exit 1
}
T1=$(date +%s)
echo "Solve+install time: $((T1-T0))s"

echo "=== Post-solve pip installs (no-deps where pkg pins old torch) ==="
# scgpt: --no-deps so it doesn't downgrade torch (it pins torch<2)
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --no-deps "scgpt>=0.2.0" || \
  echo "WARN: scgpt install failed"

# Geneformer: not on PyPI; clone first (avoiding pip's --filter=blob:none which breaks
# HF's git-LFS) then install via pip from the local checkout.
GENEFORMER_DIR=/scratch/geneformer_clone_$$
GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 \
  https://huggingface.co/ctheodoris/Geneformer "$GENEFORMER_DIR" 2>&1 | tail -5 || \
  echo "WARN: geneformer clone failed"
if [ -d "$GENEFORMER_DIR" ]; then
  micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --no-deps \
    "$GENEFORMER_DIR" || echo "WARN: geneformer local pip install failed"
  rm -rf "$GENEFORMER_DIR"
fi

# pertpy: ok with deps (no torch override in current versions)
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install pertpy || \
  echo "WARN: pertpy install partial"

# Re-pin setuptools<81 in case anything stripped pkg_resources
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --force-reinstall "setuptools<81" || \
  echo "WARN: setuptools re-pin failed"

# Best-effort flash-attn install AFTER torch is in place.
echo "=== Best-effort flash-attn install (skip on failure; SDPA fallback supported) ==="
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --no-build-isolation "flash-attn>=2.5" 2>&1 | tail -8 || \
  echo "WARN: flash-attn install failed; scGPT can run with torch SDPA fallback."

echo "=== Verification ==="
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -c "
import sys, importlib
print('Python:', sys.version)
mods = ['torch','scanpy','anndata','transformers']
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
except Exception:
    pass
for opt in ['scgpt','geneformer','flash_attn']:
    try:
        m = importlib.import_module(opt)
        print(f'  OK {opt}', getattr(m,'__version__','?'))
    except Exception as e:
        print(f'  MISSING {opt}: {e}')
"
echo "End time: $(date -Iseconds)"
echo "DONE: scGPT env build"
