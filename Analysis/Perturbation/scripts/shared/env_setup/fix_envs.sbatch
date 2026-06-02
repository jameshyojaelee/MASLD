#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=4:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=env-fix
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/env_fix_%j.log

# Fix the 4 broken/partial perturbation envs in-place via pip overrides.
# 1. perturbation_scgpt: install geneformer from git+HF (not on PyPI)
# 2. perturbation_state: replace torch with cu124 build (was cu130)
# 3. perturbation_scgen_cpa: install setuptools, replace torch, pin numpy<2, reinstall rdkit/commot
# 4. perturbation_ccc: pin numpy<2, reinstall commot

set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

ENVS_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs

eval "$(micromamba shell hook --shell bash)"

# Per-env activation + LD_LIBRARY_PATH
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

# --------------------------------------------------------------------
# 1. perturbation_scgpt — geneformer via git+HF
# --------------------------------------------------------------------
fix_env perturbation_scgpt \
  "pip install --upgrade --no-cache-dir 'git+https://huggingface.co/ctheodoris/Geneformer'" \
  "pip install --upgrade --no-cache-dir scgpt pertpy"

# --------------------------------------------------------------------
# 2. perturbation_state — replace torch with cu124 build
# --------------------------------------------------------------------
fix_env perturbation_state \
  "pip uninstall -y torch torchvision torchaudio" \
  "pip install --no-cache-dir torch==2.5.1 torchvision==0.20.1 torchaudio==2.5.1 --index-url https://download.pytorch.org/whl/cu124" \
  "pip install --upgrade --no-cache-dir 'scvi-tools>=1.1' 'lightning>=2.4,<2.6'"

# --------------------------------------------------------------------
# 3. perturbation_scgen_cpa — setuptools + torch + numpy<2
# --------------------------------------------------------------------
fix_env perturbation_scgen_cpa \
  "pip install --upgrade --no-cache-dir 'setuptools<81' wheel" \
  "pip uninstall -y torch torchvision torchaudio" \
  "pip install --no-cache-dir torch==2.4.1 torchvision==0.19.1 --index-url https://download.pytorch.org/whl/cu124" \
  "pip install --no-cache-dir 'numpy<2' --force-reinstall" \
  "pip install --no-cache-dir --force-reinstall rdkit commot" \
  "pip install --no-cache-dir --force-reinstall scgen cpa-tools pertpy"

# --------------------------------------------------------------------
# 4. perturbation_ccc — pin numpy<2 + reinstall commot
# --------------------------------------------------------------------
fix_env perturbation_ccc \
  "pip install --no-cache-dir 'numpy<2' --force-reinstall" \
  "pip install --no-cache-dir --force-reinstall commot"

# --------------------------------------------------------------------
# Final reverify across all 6 envs
# --------------------------------------------------------------------
echo ""
echo "==================== FINAL REVERIFY ===================="
bash Analysis/Perturbation/scripts/shared/env_setup/reverify_envs.sh 2>&1
echo ""
echo "[done] env_fix job complete"
