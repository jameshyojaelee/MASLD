#!/usr/bin/env bash
#SBATCH --job-name=fix_spatial
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/fix_spatial_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/fix_spatial_%j.err

module load CUDA/12.1.1

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

# Fix LD_LIBRARY_PATH
PIP_NVIDIA="${CONDA_PREFIX}/lib/python3.12/site-packages/nvidia"
for pkg in cublas cusparse cusolver curand cufft nvjitlink cusparselt; do
    d="${PIP_NVIDIA}/${pkg}/lib"
    [ -d "$d" ] && export LD_LIBRARY_PATH="${d}:${LD_LIBRARY_PATH:-}"
done

echo "=== Step 1: Install hotspot from GitHub ==="
pip install git+https://github.com/YosefLab/Hotspot.git 2>&1 || echo "WARN: Hotspot install failed"

echo ""
echo "=== Step 2: Fix RAPIDS dask/scikit-image ==="
pip install "dask==2024.12.1" "distributed==2024.12.1" "scikit-image>=0.24,<0.26" 2>&1

echo ""
echo "=== Step 3: Verify all packages ==="
python << 'PYEOF'
import sys

packages = {
    "squidpy": "squidpy",
    "cell2location": "cell2location",
    "liana": "liana",
    "hotspot": "hotspot",
    "rapids_singlecell": "rapids_singlecell",
    "dask": "dask",
    "scvi-tools": "scvi",
    "scanpy": "scanpy",
}

all_ok = True
for name, import_name in packages.items():
    try:
        mod = __import__(import_name)
        ver = getattr(mod, "__version__", "OK")
        print(f"  {name}: {ver}")
    except Exception as e:
        err_str = str(e)[:80]
        print(f"  {name}: FAILED ({err_str})")
        all_ok = False

if all_ok:
    print("\nAll packages verified successfully!")
else:
    print("\nSome packages had issues — see above")
    sys.exit(1)
PYEOF

echo ""
echo "=== Complete at $(date) ==="
