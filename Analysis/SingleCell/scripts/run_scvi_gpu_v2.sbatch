#!/usr/bin/env bash
#SBATCH --job-name=scvi_gpu_v2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/scvi_gpu_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/scvi_gpu_v2_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CONDA_ENV="${PROJECT_ROOT}/.agent/scanpy_env"
SCRIPT="${PROJECT_ROOT}/Analysis/SingleCell/scripts/scvi_integration.py"
CHECKPOINT="${PROJECT_ROOT}/Analysis/SingleCell/results_gpu/pre_scvi_checkpoint.h5ad"
OUTPUT_DIR="${PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2"

mkdir -p "${PROJECT_ROOT}/Analysis/SingleCell/logs"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/pseudobulk"
mkdir -p "${OUTPUT_DIR}/figures"

# Load CUDA 12.1 system libs (required by torch 2.2.2+cu121 for libcublas, libcudart, libcudnn, libnccl)
module load CUDA/12.1.1
module load cuDNN/8.9.2.26-CUDA-12.1.1
module load NCCL/2.18.3-GCCcore-12.3.0-CUDA-12.1.1

# CRITICAL: Module loads can set PYTHONPATH which overrides the conda env.
# Unsetting it forces Python to use the environment's site-packages first.
unset PYTHONPATH

eval "$(micromamba shell hook -s bash)"
micromamba activate "${CONDA_ENV}"

# Isolate conda env from user ~/.local packages
export PYTHONNOUSERSITE=1

# Ensure conda env libs take precedence, CUDA libs visible
export LD_LIBRARY_PATH="${CONDA_ENV}/lib:${LD_LIBRARY_PATH:-}"

# Verify environment isolation and GPU visibility
python3 -c "import numpy; print('Numpy:', numpy.__version__, '| Location:', numpy.__file__)"
python3 -c "import scipy; print('SciPy:', scipy.__version__, '| Location:', scipy.__file__)"
python3 -c "import scanpy; print('Scanpy:', scanpy.__version__, '| Location:', scanpy.__file__)"
python3 -c "import torch; print('PyTorch:', torch.__version__, '| CUDA available:', torch.cuda.is_available(), '| Device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"

# Pre-download CellTypist models if not cached
python3 -c "import celltypist; celltypist.models.download_models(force_update=False)" 2>/dev/null || true

# Unbuffered output so SLURM logs update in real-time
export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

# Run from checkpoint (skip preprocessing) into isolated output directory
python3 -X faulthandler "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --output-dir "${OUTPUT_DIR}" \
    --import-data "${CHECKPOINT}"
