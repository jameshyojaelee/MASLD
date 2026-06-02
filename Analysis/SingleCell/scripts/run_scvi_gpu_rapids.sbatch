#!/usr/bin/env bash
#SBATCH --job-name=scvi_rapids
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/scvi_rapids_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/scvi_rapids_%j.err

# scVI integration with rapids_singlecell GPU acceleration
# Uses rapids_singlecell for preprocessing + clustering (16-60x speedup)
# Uses scVI+PyTorch for model training (separate GPU acceleration)
#
# Usage:
#   sbatch run_scvi_gpu_rapids.sbatch                    # Full run from raw data
#   sbatch run_scvi_gpu_rapids.sbatch --import-data X    # Resume from checkpoint
#   sbatch run_scvi_gpu_rapids.sbatch --resume-model X   # Resume from trained model

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/Analysis/SingleCell/scripts/scvi_integration.py"
OUTPUT_DIR="${PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2"

mkdir -p "${PROJECT_ROOT}/Analysis/SingleCell/logs"
mkdir -p "${OUTPUT_DIR}"

echo "=== scVI + rapids_singlecell GPU Pipeline ==="
echo "Node: $(hostname)"
echo "Date: $(date)"
nvidia-smi
echo ""

# Activate rapids_singlecell environment (has RAPIDS + scanpy + scvi)
eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

# Add pip-installed NVIDIA libs to LD_LIBRARY_PATH (cusparselt for torch)
SITE_PKGS="${CONDA_PREFIX}/lib/python3.12/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

# Verify GPU
python3 -c "
import cupy as cp
print(f'CuPy: {cp.__version__}, CUDA: {cp.cuda.is_available()}, GPU: {cp.cuda.runtime.getDeviceCount()}')
import rapids_singlecell as rsc
print(f'rapids_singlecell: {rsc.__version__}')
"

# Run with GPU acceleration (--use-gpu enables rapids_singlecell)
python3 -X faulthandler "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --output-dir "${OUTPUT_DIR}" \
    --use-gpu \
    "$@"

echo ""
echo "Done: $(date)"
