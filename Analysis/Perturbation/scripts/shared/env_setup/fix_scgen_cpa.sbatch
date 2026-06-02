#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=4:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=scGEN
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_scgen_cpa_rebuild_%j.log

# Rebuild perturbation_scgen_cpa from scratch — torch 1.13 keeps being re-installed
# by scvi-tools/cpa pip deps, undoing the cu124 pin. Clean wipe + channel-namespaced
# specs + explicit constraints on torch in pip install lines.

set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_scgen_cpa

echo "[+] Wiping broken env at $ENV_PATH"
rm -rf "$ENV_PATH"

echo "[+] Building env from YAML"
eval "$(micromamba shell hook --shell bash)"
micromamba env create -f envs/perturbation/perturbation_scgen_cpa.yml -p "$ENV_PATH" -y

# Post-solve: use pip with strict torch constraint to prevent downgrade
export LD_LIBRARY_PATH="$ENV_PATH/lib:${LD_LIBRARY_PATH:-}"
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 bash -c '
  pip install --upgrade pip "setuptools<81" wheel
  # Force torch from pytorch index, lock to 2.4.1+cu124
  pip install --no-cache-dir "torch==2.4.1" "torchvision==0.19.1" --index-url https://download.pytorch.org/whl/cu124
  # Pin numpy<2 BEFORE installing pkgs that pull it in
  pip install --no-cache-dir --upgrade "numpy<2"
  # Now install scgen / cpa / pertpy with --no-deps so they cannot drag in different torch
  pip install --no-cache-dir --no-deps scgen
  pip install --no-cache-dir --no-deps cpa-tools
  pip install --no-cache-dir --no-deps pertpy
  # Reinstall their non-torch deps explicitly (the ones we want)
  pip install --no-cache-dir --no-deps mudata rdkit "numpy<2"
'

# Verify
echo "[+] Smoke test"
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -c "
import torch; print('torch:', torch.__version__, 'cuda:', torch.version.cuda)
import numpy as np; print('numpy:', np.__version__)
import scvi; print('scvi:', scvi.__version__)
import scgen; print('scgen:', scgen.__version__)
try:
    import cpa; print('cpa: imported OK')
except Exception as e: print('cpa FAIL:', str(e)[:100])
try:
    from rdkit import Chem; print('rdkit: imported OK')
except Exception as e: print('rdkit FAIL:', str(e)[:100])
"
echo "[done] scGEN+CPA env rebuild complete"
