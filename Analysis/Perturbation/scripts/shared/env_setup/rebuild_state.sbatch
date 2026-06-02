#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=4:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=STATE
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_state_rebuild_%j.log

# Rebuild perturbation_state from scratch — previous env destroyed by env-fix.
# Channel-namespaced specs to force pytorch channel to win over conda-forge.

set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

ENV_PATH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_state

# Wipe broken env
echo "[+] Removing broken env at $ENV_PATH"
rm -rf "$ENV_PATH"

# Build from YAML
echo "[+] Building env from envs/perturbation/perturbation_state.yml"
eval "$(micromamba shell hook --shell bash)"
micromamba env create -f envs/perturbation/perturbation_state.yml -p "$ENV_PATH" -y

# Post-solve: install STATE via uv tool
echo "[+] Installing arc-state via pip (uv unavailable, use pip)"
export LD_LIBRARY_PATH="$ENV_PATH/lib:${LD_LIBRARY_PATH:-}"
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 bash -c '
  pip install --upgrade pip
  pip install --no-cache-dir arc-state || pip install --no-cache-dir "git+https://github.com/ArcInstitute/state.git" || echo "STATE install via pip failed — install manually later"
'

# Verify
echo "[+] Smoke test"
micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -c "
import torch
print('torch:', torch.__version__, 'cuda:', torch.version.cuda)
import scvi; print('scvi:', scvi.__version__)
import lightning; print('lightning:', lightning.__version__)
import transformers; print('transformers:', transformers.__version__)
"
echo "[done] STATE env rebuild complete"
