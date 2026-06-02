#!/usr/bin/env bash
#SBATCH --job-name=test_gpu_pipeline
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/test_gpu/test_gpu_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/test_gpu/test_gpu_%j.err

# End-to-end GPU pipeline test
# Tests: per-sample scanpy, multi-sample preprocessing, mini scVI
# All output in test_gpu/ — no production files touched

set -euo pipefail

TEST_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/test_gpu"

echo "=== GPU Pipeline Test ==="
echo "Node: $(hostname) | Date: $(date)"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

# Add pip-installed NVIDIA libs to LD_LIBRARY_PATH (cusparselt for torch)
SITE_PKGS="${CONDA_PREFIX}/lib/python3.12/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

# Verify imports
python3 -c "
import cupy as cp
print(f'CuPy: {cp.__version__}, GPU: {cp.cuda.runtime.getDeviceCount()}')
import rapids_singlecell as rsc
print(f'rapids_singlecell: {rsc.__version__}')
import scvi
print(f'scvi-tools: {scvi.__version__}')
import scanpy as sc
print(f'scanpy: {sc.__version__}')
"

echo ""
echo "--- Running tests (GPU + CPU comparison) ---"
echo ""

python3 -X faulthandler "${TEST_DIR}/test_gpu_pipeline.py" \
    --output-dir "${TEST_DIR}" \
    --tests 1 2 3

echo ""
echo "Done: $(date)"
