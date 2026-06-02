#!/usr/bin/env bash
#SBATCH --job-name=test_scalesc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/test_scalesc/test_scalesc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/test_scalesc/test_scalesc_%j.err

# ScaleSC GPU pipeline test
# Tests: full pipeline on Liver Atlas samples + comparison vs rapids_singlecell

set -euo pipefail

TEST_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/test_scalesc"

echo "=== ScaleSC GPU Pipeline Test ==="
echo "Node: $(hostname) | Date: $(date)"
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv,noheader
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

export PYTHONUNBUFFERED=1
export PYTHONFAULTHANDLER=1

# Add pip-installed NVIDIA libs to LD_LIBRARY_PATH
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
import scalesc as ssc
print(f'ScaleSC: {ssc.__version__}')
import scanpy as sc
print(f'scanpy: {sc.__version__}')
"

echo ""
echo "--- Running ScaleSC test (5 samples) ---"
echo ""

# Test with 5 samples first (manageable size, ~50K cells)
python3 -X faulthandler "${TEST_DIR}/test_scalesc.py" \
    --output-dir "${TEST_DIR}/output_5samples" \
    --n-samples 5 \
    --max-cell-batch 100000

echo ""
echo "--- Running ScaleSC test (all 48 samples) ---"
echo ""

# Full test with all 48 samples (~500K cells)
python3 -X faulthandler "${TEST_DIR}/test_scalesc.py" \
    --output-dir "${TEST_DIR}/output_full" \
    --max-cell-batch 100000 \
    --skip-baseline

echo ""
echo "Done: $(date)"
