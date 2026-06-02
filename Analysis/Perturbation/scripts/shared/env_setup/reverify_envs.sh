#!/bin/bash
# reverify_envs.sh — standalone verification for the 6 perturbation envs.
# Usage: bash Analysis/Perturbation/scripts/shared/env_setup/reverify_envs.sh
# Or run a single env:  bash reverify_envs.sh perturbation_state
#
# CRITICAL: this script sets PYTHONNOUSERSITE=1 + LD_LIBRARY_PATH=$CONDA_PREFIX/lib
# so that user-site ~/.local packages do NOT leak in. If you forget this, you will
# see "OK torch 2.1.2+cu121" which is the user-site leak, not the env's pytorch.

set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

ENVS_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs

REQUESTED="${1:-all}"
if [ "$REQUESTED" = "all" ]; then
  ENVS=(perturbation_gears perturbation_scgen_cpa perturbation_state perturbation_tahoe perturbation_scgpt perturbation_ccc)
else
  ENVS=("$REQUESTED")
fi

# Per-env import smoke tests. Add OK/FAIL marker on each.
declare -A SMOKE
SMOKE[perturbation_gears]="torch scanpy anndata networkx pertpy llvmlite gears"
SMOKE[perturbation_scgen_cpa]="torch scvi lightning pertpy scgen cpa ot"
SMOKE[perturbation_state]="torch scvi lightning transformers pertpy"
SMOKE[perturbation_tahoe]="torch transformers huggingface_hub rdkit pertpy"
SMOKE[perturbation_scgpt]="torch scanpy transformers scgpt geneformer pertpy"
SMOKE[perturbation_ccc]="torch scanpy squidpy omnipath ot liana commot pertpy"

for E in "${ENVS[@]}"; do
  ENV_PATH="$ENVS_DIR/$E"
  echo "==================== $E ===================="
  if [ ! -d "$ENV_PATH" ]; then
    echo "  MISSING env at $ENV_PATH"
    continue
  fi
  export CONDA_PREFIX="$ENV_PATH"
  export LD_LIBRARY_PATH="$ENV_PATH/lib:${LD_LIBRARY_PATH:-}"
  MODS="${SMOKE[$E]:-torch scanpy anndata}"
  micromamba run -p "$ENV_PATH" \
    env PYTHONNOUSERSITE=1 \
    python -c "
import sys, importlib
print('Python exe :', sys.executable)
print('Python ver :', sys.version.split()[0])
print('Prefix     :', sys.prefix)
print('User-site  : ENABLED' if (hasattr(sys,'flags') and not sys.flags.no_user_site)==False else 'User-site  : DISABLED (good)')
mods = '''$MODS'''.split()
n_ok=0; n_fail=0
for m in mods:
    try:
        mod = importlib.import_module(m)
        v = getattr(mod, '__version__', '?')
        loc = getattr(mod,'__file__','?')
        ok_path = loc.startswith('$ENV_PATH')
        print(f'  {\"OK \" if ok_path else \"OK*\"} {m} {v}  ({loc[:80]})')
        n_ok+=1
    except Exception as e:
        print(f'  FAIL {m}: {type(e).__name__}: {str(e)[:160]}')
        n_fail+=1
try:
    import torch
    print('  CUDA available:', torch.cuda.is_available())
    print('  torch.version.cuda:', getattr(torch.version,'cuda',None))
except Exception as e:
    print('  torch CUDA check failed:', e)
print(f'  SUMMARY: ok={n_ok} fail={n_fail}')
"
done
echo "==================== DONE ===================="
