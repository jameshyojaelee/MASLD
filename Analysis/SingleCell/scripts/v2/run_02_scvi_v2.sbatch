#!/usr/bin/env bash
#SBATCH --job-name=scvi_v2
#SBATCH --partition=gpu
#SBATCH --qos=interactive
#SBATCH --gres=gpu:1
#SBATCH --exclude=ne1dg7-[001-010]
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2_phase05/logs/scvi_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2_phase05/logs/scvi_v2_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CONDA_ENV="${PROJECT_ROOT}/.agent/scanpy_env"
SCRIPT="${PROJECT_ROOT}/Analysis/SingleCell/scripts/v2/02_train_scvi_v2.py"

mkdir -p "${PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2_phase05/logs"

# Load CUDA system libs to match torch 2.3.1+cu121
module load CUDA/12.1.1 || true
module load cuDNN/8.9.2.26-CUDA-12.1.1 || true

# Clear PYTHONPATH to avoid GCCcore libstdc++ vs env libs mismatch
unset PYTHONPATH
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="${CONDA_ENV}/lib:${LD_LIBRARY_PATH:-}"

# scvi-tools network / telemetry off
export WANDB_MODE=disabled
export WANDB_DISABLED=true
export SCVI_NO_TELEMETRY=true

export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

# Diagnostics
"${CONDA_ENV}/bin/python" -c "
import torch, scvi, anndata, scanpy
print('scvi-tools:', scvi.__version__)
print('anndata:   ', anndata.__version__)
print('scanpy:    ', scanpy.__version__)
print('torch:     ', torch.__version__)
print('cuda:      ', torch.cuda.is_available())
if torch.cuda.is_available():
    print('device:    ', torch.cuda.get_device_name(0))
"

"${CONDA_ENV}/bin/python" -X faulthandler "${SCRIPT}"
