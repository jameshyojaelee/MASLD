#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=2:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=env-fix-3
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_fix3_%j.log

# Targeted fixes for the 3 envs that still have missing/broken packages
# after env-doctor's rebuild. The torch + numpy + nccl base is correct;
# just need to add missing high-level pkgs.
#
# scGEN env: pertpy, rdkit, scgen-2.x (scvi 1.4 API compat)
# scGPT env: IPython, geneformer (git+HF)
# CCC env: zarr<3 (for squidpy), reinstall commot, reinstall squidpy

set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

ENVS_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs
eval "$(micromamba shell hook --shell bash)"

fix_env() {
  local ENV_NAME=$1
  shift
  export CONDA_PREFIX="$ENVS_DIR/$ENV_NAME"
  export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
  echo "==================== Fixing $ENV_NAME ===================="
  for cmd in "$@"; do
    echo "[+] $cmd"
    micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 bash -c "$cmd" || echo "  (continuing despite error)"
  done
  echo ""
}

# ----------------------------------------------------------------------
# scGEN env: install pertpy + rdkit via conda, scgen via pip --no-deps
# ----------------------------------------------------------------------
fix_env perturbation_scgen_cpa \
  "micromamba install -p $ENVS_DIR/perturbation_scgen_cpa -c conda-forge -y rdkit pertpy" \
  "pip install --no-cache-dir --no-deps 'scgen>=2.2,<3' || pip install --no-cache-dir --no-deps scgen" \
  "pip install --no-cache-dir --no-deps cpa-tools"

# ----------------------------------------------------------------------
# scGPT env: IPython via conda, Geneformer via git+HF
# ----------------------------------------------------------------------
fix_env perturbation_scgpt \
  "micromamba install -p $ENVS_DIR/perturbation_scgpt -c conda-forge -y ipython" \
  "pip install --no-cache-dir 'git+https://huggingface.co/ctheodoris/Geneformer'"

# ----------------------------------------------------------------------
# CCC env: pin zarr<3 (for squidpy compat) + reinstall squidpy + commot
# ----------------------------------------------------------------------
fix_env perturbation_ccc \
  "pip install --no-cache-dir --force-reinstall 'zarr<3' 'numpy<2'" \
  "pip install --no-cache-dir --force-reinstall --no-deps squidpy" \
  "pip install --no-cache-dir --force-reinstall --no-deps commot"

# ----------------------------------------------------------------------
# Final verify
# ----------------------------------------------------------------------
echo ""
echo "==================== FINAL REVERIFY (3 fixed envs) ===================="
for E in perturbation_scgen_cpa perturbation_scgpt perturbation_ccc; do
  export CONDA_PREFIX="$ENVS_DIR/$E"
  export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
  echo "--- $E ---"
  micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 python -c "
import importlib
mods = {
  'perturbation_scgen_cpa': 'torch scanpy anndata scvi lightning pertpy scgen cpa rdkit ot',
  'perturbation_scgpt':     'torch scanpy transformers scgpt geneformer pertpy',
  'perturbation_ccc':       'torch scanpy squidpy omnipath ot liana commot pertpy',
}['$E'].split()
ok, fail = 0, 0
for m in mods:
  try:
    mod = importlib.import_module(m)
    v = getattr(mod, '__version__', '?')
    print(f'  OK  {m} {v}'); ok += 1
  except Exception as e:
    print(f'  FAIL {m}: {type(e).__name__}: {str(e)[:120]}'); fail += 1
print(f'  SUMMARY: ok={ok} fail={fail}')
"
  echo ""
done
echo "[done] env_fix_3 complete"
