#!/bin/bash
#SBATCH --job-name=geneformer-inst
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/install_geneformer_%j.log

# Targeted Geneformer install in scGPT env (git+HF approach was failing in
# env-final due to git LFS + auth in batch context). Use huggingface_hub to
# clone via Python API which handles auth + LFS more gracefully.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1

eval "$(micromamba shell hook --shell bash)"
micromamba activate perturbation_scgpt

set -u

# Use the Geneformer checkpoint we already downloaded at data/perturbation/checkpoints/geneformer/Geneformer/
# Install the geneformer package via pip from that local dir
GENEFORMER_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/checkpoints/geneformer/Geneformer
if [ -d "$GENEFORMER_DIR" ]; then
  echo "[+] Installing geneformer from local checkpoint dir $GENEFORMER_DIR"
  ls "$GENEFORMER_DIR"/ | head
  # The HF repo IS a python package: setup.py / pyproject.toml at root
  if [ -f "$GENEFORMER_DIR/setup.py" ] || [ -f "$GENEFORMER_DIR/pyproject.toml" ]; then
    pip install --no-cache-dir "$GENEFORMER_DIR"
  else
    echo "[!] No setup.py/pyproject.toml in $GENEFORMER_DIR — installing via git+HF"
    pip install --no-cache-dir 'git+https://huggingface.co/ctheodoris/Geneformer'
  fi
else
  echo "[!] Local Geneformer dir missing; falling back to git+HF"
  pip install --no-cache-dir 'git+https://huggingface.co/ctheodoris/Geneformer'
fi

# Also install IPython (case-insensitive — verify script tries lowercase `ipython`)
pip install --no-cache-dir ipython || echo "ipython already present"

# Pertpy (verify failed last round — try via conda since pip is fragile)
micromamba install -p .envs/perturbation_scgpt -c conda-forge -y pertpy || pip install --no-cache-dir pertpy

# Pydantic (needed for runner_base output_schema validation)
micromamba install -p .envs/perturbation_scgpt -c conda-forge -y pydantic || pip install --no-cache-dir pydantic

echo ""
echo "[verify]"
python -c "
import importlib
for m in 'torch scanpy transformers scgpt geneformer pertpy IPython pydantic'.split():
    try:
        mod = importlib.import_module(m)
        v = getattr(mod, '__version__', '?')
        print(f'  OK  {m} {v}')
    except Exception as e:
        print(f'  FAIL {m}: {type(e).__name__}: {str(e)[:80]}')
"

echo "[done] install-geneformer complete"
