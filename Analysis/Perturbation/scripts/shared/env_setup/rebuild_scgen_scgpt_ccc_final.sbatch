#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=6:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=env-final
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_final_%j.log

# Final attempt: rebuild scGEN + scGPT from scratch sequentially (avoid mamba
# lock contention), then surgical fix for CCC numpy conflict.
#
# Why sequential: the previous parallel rebuild + parallel pip fix caused env
# destruction (scGEN/scGPT lost their python binary mid-stream when other
# install jobs touched the shared pkgs cache).

set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

ENVS_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs
eval "$(micromamba shell hook --shell bash)"

# ----------------------------------------------------------------------
# 1. scGEN — full rebuild from YAML, then post-install via pip (NO --no-deps
#    chaos this time; let pip resolve dependencies properly).
# ----------------------------------------------------------------------
echo "==================== scGEN rebuild ===================="
rm -rf "$ENVS_DIR/perturbation_scgen_cpa"
micromamba env create -f envs/perturbation/perturbation_scgen_cpa.yml -p "$ENVS_DIR/perturbation_scgen_cpa" -y

export CONDA_PREFIX="$ENVS_DIR/perturbation_scgen_cpa"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

# Verify base then install scgen via pip (scgen needs scvi-tools, which conda installed)
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 python -c "
import torch; print('torch:', torch.__version__, 'cuda:', torch.version.cuda)
import scvi; print('scvi:', scvi.__version__)
"
# scgen 2.1.x requires old scvi._compat API; latest scvi removed it. Use scvi.external.SCGEN instead.
# pertpy + rdkit installed via conda in env-fix-3 but env was destroyed; reinstall.
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 bash -c '
  micromamba install -p '"$CONDA_PREFIX"' -c conda-forge -y rdkit pertpy
'
# scvi.external.SCGEN is bundled with scvi-tools >= 1.0 — no separate scgen package needed
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 python -c "
from scvi.external import SCGEN; print('SCGEN (from scvi.external):', SCGEN)
"

# CPA via pip
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 pip install --no-cache-dir cpa-tools || echo "cpa-tools install failed"

echo "[scGEN] verify ..."
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 python -c "
import importlib
for m in 'torch scanpy anndata scvi lightning pertpy rdkit ot cpa'.split():
    try:
        mod = importlib.import_module(m)
        print(f'  OK  {m} {getattr(mod, \"__version__\", \"?\")}')
    except Exception as e: print(f'  FAIL {m}: {type(e).__name__}: {str(e)[:80]}')
try:
    from scvi.external import SCGEN; print('  OK  scvi.external.SCGEN')
except Exception as e: print(f'  FAIL SCGEN: {e}')
"

# ----------------------------------------------------------------------
# 2. scGPT — full rebuild from YAML, then post-install scGPT + Geneformer + IPython
# ----------------------------------------------------------------------
echo ""
echo "==================== scGPT rebuild ===================="
rm -rf "$ENVS_DIR/perturbation_scgpt"
micromamba env create -f envs/perturbation/perturbation_scgpt.yml -p "$ENVS_DIR/perturbation_scgpt" -y

export CONDA_PREFIX="$ENVS_DIR/perturbation_scgpt"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 python -c "
import torch; print('torch:', torch.__version__, 'cuda:', torch.version.cuda)
"
# IPython via conda
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 bash -c '
  micromamba install -p '"$CONDA_PREFIX"' -c conda-forge -y ipython
'
# scGPT via pip
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 pip install --no-cache-dir scgpt || echo "scgpt install failed"
# Geneformer via git clone + pip (NOT git+https — sometimes HF clone is more reliable)
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 bash -c '
  cd /tmp
  rm -rf Geneformer_install
  git clone https://huggingface.co/ctheodoris/Geneformer Geneformer_install || echo "Geneformer clone failed"
  if [ -d Geneformer_install ]; then
    cd Geneformer_install
    pip install --no-cache-dir . || echo "Geneformer install failed"
  fi
'
# pertpy
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 bash -c '
  micromamba install -p '"$CONDA_PREFIX"' -c conda-forge -y pertpy
'

echo "[scGPT] verify ..."
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 python -c "
import importlib
for m in 'torch scanpy transformers scgpt geneformer pertpy ipython'.split():
    try:
        mod = importlib.import_module(m)
        print(f'  OK  {m} {getattr(mod, \"__version__\", \"?\")}')
    except Exception as e: print(f'  FAIL {m}: {type(e).__name__}: {str(e)[:80]}')
"

# ----------------------------------------------------------------------
# 3. CCC — surgical fix for numpy conflict
#    commot needs np.Inf (numpy<2)
#    ot/pertpy/squidpy need StringDType (numpy>=2)
#    Resolution: pin numpy 1.26.4 (last 1.x release w/ broad compatibility);
#    upgrade commot if available (commot >=0.1 supports numpy 2); else accept
#    commot as broken and use NICHES alone for D4.
# ----------------------------------------------------------------------
echo ""
echo "==================== CCC surgical fix ===================="
export CONDA_PREFIX="$ENVS_DIR/perturbation_ccc"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

# Try upgrading commot to a version that supports numpy 2
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 bash -c '
  pip install --no-cache-dir --upgrade "commot>=0.1" || echo "no commot>=0.1 — staying on old version"
'
# Force numpy to a working middle ground (1.26.4 is widely compatible)
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 pip install --no-cache-dir --force-reinstall "numpy==1.26.4"
# Reinstall ot/pertpy with numpy 1.26 compat
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 pip install --no-cache-dir --force-reinstall --no-deps "pot<0.10" pertpy

echo "[CCC] verify ..."
micromamba run -p "$CONDA_PREFIX" env PYTHONNOUSERSITE=1 python -c "
import importlib
for m in 'torch scanpy squidpy omnipath ot liana commot pertpy'.split():
    try:
        mod = importlib.import_module(m)
        print(f'  OK  {m} {getattr(mod, \"__version__\", \"?\")}')
    except Exception as e: print(f'  FAIL {m}: {type(e).__name__}: {str(e)[:80]}')
"

echo ""
echo "[done] env_final complete"
