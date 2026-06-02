#!/usr/bin/env bash
#SBATCH --job-name=install_spatial
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/install_spatial_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/install_spatial_%j.err
##############################################################################
# Install spatial analysis packages into rapids_singlecell environment.
#
# MUST run on GPU node because squidpy→spatialdata→datashader→cudf
# import chain requires CUDA at import time.
#
# After install: fix RAPIDS-broken dask/scikit-image versions.
##############################################################################
set -euo pipefail

module load CUDA/12.1.1

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

# Fix LD_LIBRARY_PATH for pip-installed NVIDIA libs
PIP_NVIDIA="${CONDA_PREFIX}/lib/python3.12/site-packages/nvidia"
for pkg in cublas cusparse cusolver curand cufft nvjitlink cusparselt; do
    d="${PIP_NVIDIA}/${pkg}/lib"
    [ -d "$d" ] && export LD_LIBRARY_PATH="${d}:${LD_LIBRARY_PATH:-}"
done

echo "=== Installing spatial packages into rapids_singlecell ==="
echo "Python: $(which python)"
echo "Env: ${CONDA_PREFIX}"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"

# Step 1: Install all packages (don't verify between installs to avoid partial failures)
echo ""
echo "--- Installing cell2location ---"
pip install cell2location 2>&1

echo ""
echo "--- Installing liana-py ---"
pip install liana 2>&1

echo ""
echo "--- Installing hotspot ---"
pip install hotspot 2>&1

# Step 2: Fix RAPIDS dask/scikit-image breakage
# squidpy 1.8.1 pulled dask 2026.x and scikit-image 0.26.0
# RAPIDS 25.02 requires dask==2024.12.1 and scikit-image<0.26
echo ""
echo "--- Fixing RAPIDS-compatible versions ---"
pip install "dask==2024.12.1" "distributed==2024.12.1" "scikit-image>=0.24,<0.26" 2>&1

# Step 3: Verify all packages (on GPU node, so cudf import should work)
echo ""
echo "=== Verification ==="
python << 'PYEOF'
results = {}
for pkg_name in ["squidpy", "cell2location", "liana", "hotspot"]:
    try:
        mod = __import__(pkg_name)
        ver = getattr(mod, "__version__", "OK")
        results[pkg_name] = f"{ver}"
        print(f"  {pkg_name}: {ver}")
    except Exception as e:
        results[pkg_name] = f"FAILED: {e}"
        print(f"  {pkg_name}: FAILED ({e})")

# Verify RAPIDS still works
try:
    import rapids_singlecell
    print(f"  rapids_singlecell: {rapids_singlecell.__version__}")
except Exception as e:
    print(f"  rapids_singlecell: FAILED ({e})")

try:
    import dask
    print(f"  dask: {dask.__version__}")
except Exception as e:
    print(f"  dask: FAILED ({e})")

try:
    import scvi
    print(f"  scvi-tools: {scvi.__version__}")
except Exception as e:
    print(f"  scvi-tools: FAILED ({e})")

all_ok = all("FAILED" not in v for v in results.values())
print(f"\n{'All spatial packages OK!' if all_ok else 'Some packages failed — see above'}")
PYEOF

echo ""
echo "=== Complete at $(date) ==="
