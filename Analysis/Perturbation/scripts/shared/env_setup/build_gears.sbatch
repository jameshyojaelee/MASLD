#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=env_gears
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_gears.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Prevent user-site (~/.local/lib/pythonX.Y/site-packages) from leaking into env Python.
# Without this, stale torch 2.1.2+cu121 and numpy 2.2.6 from ~/.local shadow conda packages.
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_gears
YAML=envs/perturbation/perturbation_gears.yml
# LD_LIBRARY_PATH set after env exists (uses ENV_PATH as CONDA_PREFIX surrogate)
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
  echo "ERROR: micromamba env create failed for GEARS env"
  exit 1
}
T1=$(date +%s)
echo "Solve+install time: $((T1-T0))s"

echo "=== Verification ==="
micromamba run -p "$ENV_PATH" python -c "
import sys, importlib
print('Python:', sys.version)
mods = ['torch','torch_geometric','scanpy','anndata','networkx','pertpy']
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
try:
    import gears
    print('  OK gears', getattr(gears,'__version__','?'))
except Exception as e:
    print(f'  GEARS NOT IMPORTABLE: {e}')
"
echo "End time: $(date -Iseconds)"
echo "DONE: GEARS env build"
