#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=env_scgen_cpa
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_scgen_cpa.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Prevent user-site (~/.local/lib/pythonX.Y/site-packages) leak into env Python.
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_scgen_cpa
YAML=envs/perturbation/perturbation_scgen_cpa.yml
export LD_LIBRARY_PATH="$ENV_PATH/lib:${LD_LIBRARY_PATH:-}"

echo "=========================================="
echo "Building $ENV_PATH from $YAML"
echo "Start time: $(date -Iseconds)"
echo "=========================================="

if [ -d "$ENV_PATH" ]; then
  echo "Removing existing partial env at $ENV_PATH"
  rm -rf "$ENV_PATH"
fi

T0=$(date +%s)
micromamba env create -f "$YAML" -p "$ENV_PATH" -y || {
  echo "ERROR: micromamba env create failed for scGen/CPA env"
  exit 1
}
T1=$(date +%s)
echo "Solve+install time: $((T1-T0))s"

echo "=== Post-solve pip install (with --no-deps where needed) ==="
# CRITICAL: --no-deps so cpa-tools doesn't downgrade torch to 1.13.1+cu117
#           and scgen doesn't pin torch==2.1.2, replacing conda's 2.4.1.
# NOTE: standalone scgen 2.1.0 imports scvi._compat which was removed in scvi-tools 1.x.
# Strategy: SKIP standalone scgen — use scvi.external.SCGEN from scvi-tools 1.4 directly.
# Only install cpa-tools + pertpy + ott-jax with --no-deps.
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --no-deps \
  "cpa-tools>=0.7" ott-jax || \
  echo "WARN: cpa-tools install partial"

# pertpy: install with deps but block torch/numpy/setuptools upgrades via constraints.
# Use a real tmp file (not <() process substitution — doesn't work in SLURM dash shell).
CONSTRAINTS_FILE=$(mktemp /tmp/scgen_cpa_constraints.XXXXXX.txt)
cat > "$CONSTRAINTS_FILE" <<EOF
torch==2.4.1
numpy<2
setuptools<81
lightning<2.6
EOF
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install \
  --constraint "$CONSTRAINTS_FILE" \
  pertpy || \
  echo "WARN: pertpy install partial (constraint may not have prevented all overrides)"
rm -f "$CONSTRAINTS_FILE"

# Re-install setuptools<81 to be safe
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -m pip install --force-reinstall \
  "setuptools<81" || echo "WARN: setuptools re-pin failed"

echo "=== Verification ==="
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -c "
import sys, importlib
print('Python:', sys.version)
mods = ['torch','scanpy','anndata','scvi','lightning','pertpy']
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
for opt in ['cpa','ot']:
    try:
        m = importlib.import_module(opt)
        print(f'  OK {opt}', getattr(m,'__version__','?'))
    except Exception as e:
        print(f'  MISSING {opt}: {e}')
# scGen via scvi-tools 1.4 internal model (standalone scgen 2.1 broken on modern scvi)
try:
    from scvi.external import SCGEN
    print('  OK scvi.external.SCGEN (use this instead of standalone scgen pkg)')
except Exception as e:
    print(f'  MISSING scvi.external.SCGEN: {e}')
"
echo "End time: $(date -Iseconds)"
echo "DONE: scGen+CPA env build"
