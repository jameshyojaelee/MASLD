#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=8:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=env_tahoe
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_tahoe.log

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Prevent user-site (~/.local/lib/pythonX.Y/site-packages) leak into env Python.
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_tahoe
YAML=envs/perturbation/perturbation_tahoe.yml
export LD_LIBRARY_PATH="$ENV_PATH/lib:${LD_LIBRARY_PATH:-}"

echo "=========================================="
echo "Building $ENV_PATH from $YAML"
echo "Start time: $(date -Iseconds)"
echo "NOTE: vevo-sdk PyPI placeholder may not exist; install will be tried with fallbacks"
echo "=========================================="

if [ -d "$ENV_PATH" ]; then
  echo "Removing existing partial env at $ENV_PATH"
  rm -rf "$ENV_PATH"
fi

# Patch YAML to strip placeholder pip pkg that may fail to resolve;
# we install vevo SDK candidates manually after env solve.
PATCHED_YAML=envs/perturbation/perturbation_tahoe_patched.yml
python3 - <<PY
import re
src = open("$YAML").read()
# Drop the placeholder "- vevo-sdk" line so the solve succeeds; install via pip after.
src = re.sub(r"^\s*-\s*vevo-sdk.*$", "    # vevo-sdk attempted post-install", src, flags=re.M)
open("$PATCHED_YAML","w").write(src)
print("Patched YAML written:", "$PATCHED_YAML")
PY

T0=$(date +%s)
micromamba env create -f "$PATCHED_YAML" -p "$ENV_PATH" -y || {
  echo "ERROR: micromamba env create failed for Tahoe env"
  exit 1
}
T1=$(date +%s)
echo "Solve+install time: $((T1-T0))s"

echo "=== Try Vevo/Tahoe SDK pip install candidates ==="
SDK_INSTALLED=""
for pkg in vevo-state vevo-sdk vevotx tahoebio tahoe-x1 tahoe-100m; do
  if micromamba run -p "$ENV_PATH" python -m pip install "$pkg" 2>&1 | tee /tmp/tahoe_$pkg.log | grep -qE "Successfully installed"; then
    SDK_INSTALLED="$pkg"
    echo "  INSTALLED: $pkg"
    break
  else
    echo "  not on PyPI: $pkg"
  fi
done

if [ -z "$SDK_INSTALLED" ]; then
  echo "No PyPI package found. Attempting git+huggingface install for Tahoe-x1 model card repo."
  micromamba run -p "$ENV_PATH" python -m pip install \
    git+https://huggingface.co/tahoebio/Tahoe-x1 2>&1 | tail -5 || true
fi

echo "=== Verification ==="
micromamba run -p "$ENV_PATH" python -c "
import sys, importlib
print('Python:', sys.version)
mods = ['torch','scanpy','anndata','transformers','huggingface_hub','rdkit','pertpy']
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
for opt in ['vevo','vevo_sdk','vevotx','tahoebio','tahoe_x1']:
    try:
        m = importlib.import_module(opt)
        print(f'  OK {opt}', getattr(m,'__version__','?'))
    except Exception:
        pass
print('Note: Tahoe weights downloaded separately via huggingface_hub.snapshot_download from datasets/tahoebio/Tahoe-100M')
"
echo "End time: $(date -Iseconds)"
echo "DONE: Tahoe env build"
