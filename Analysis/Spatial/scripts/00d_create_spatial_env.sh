#!/usr/bin/env bash
#SBATCH --job-name=spatial_env2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/spatial_env2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/spatial_env2_%j.err
set -euo pipefail

module load CUDA/12.1.1

eval "$(micromamba shell hook -s bash)"

# Remove existing spatial env
micromamba env remove -n spatial -y 2>/dev/null || true

echo "=== Creating spatial environment ==="
micromamba create -n spatial python=3.12 -y -c conda-forge 2>&1 | tail -3

micromamba activate spatial

echo ""
echo "=== Installing PyTorch (CUDA 12.4) ==="
pip install torch==2.6.0+cu124 --index-url https://download.pytorch.org/whl/cu124 --no-deps 2>&1 | tail -3

echo ""
echo "=== Installing core packages (excluding hotspot) ==="
pip install \
    "scanpy>=1.10" \
    "anndata>=0.11" \
    "squidpy>=1.6,<2.0" \
    "cell2location>=0.1.5" \
    "liana>=1.7" \
    "scvi-tools>=1.4" \
    "spatialdata>=0.2" \
    "dask>=2025.2,<2026" \
    "distributed>=2025.2,<2026" \
    "pyro-ppl>=1.9" \
    "networkx>=3.0" \
    "matplotlib>=3.8" \
    "seaborn>=0.13" \
    "pyyaml" \
    "pooch" \
    "matplotlib-scalebar" \
    "validators" \
    "omnipath" \
    2>&1 | tail -15

echo ""
echo "=== Installing hotspot from GitHub ==="
pip install "git+https://github.com/YosefLab/Hotspot.git" 2>&1 | tail -5

echo ""
echo "=== Verify all packages ==="
python << 'PYEOF'
import traceback

results = {}
packages = [
    ("squidpy", "squidpy"),
    ("cell2location", "cell2location"),
    ("liana", "liana"),
    ("hotspot", "hotspot"),
    ("scanpy", "scanpy"),
    ("scvi-tools", "scvi"),
    ("anndata", "anndata"),
    ("dask", "dask"),
    ("spatialdata", "spatialdata"),
    ("torch", "torch"),
    ("pyro-ppl", "pyro"),
    ("networkx", "networkx"),
]

for name, imp in packages:
    try:
        mod = __import__(imp)
        ver = getattr(mod, "__version__", "OK")
        results[name] = ver
        print(f"  {name}: {ver}")
    except Exception as e:
        results[name] = "FAILED"
        print(f"  {name}: FAILED ({str(e)[:100]})")

# Test squidpy spatial functions
print("\nTesting squidpy.gr functions...")
try:
    from squidpy import gr
    for f in ["spatial_neighbors", "spatial_autocorr", "nhood_enrichment", "ligrec"]:
        assert hasattr(gr, f), f"{f} missing"
    print("  All squidpy.gr spatial analysis functions: OK")
except Exception as e:
    print(f"  squidpy.gr: FAILED ({e})")
    traceback.print_exc()

# Test cell2location model import
print("\nTesting cell2location models...")
try:
    from cell2location.models import RegressionModel, Cell2location
    print("  cell2location models: OK")
except Exception as e:
    print(f"  cell2location models: FAILED ({e})")
    traceback.print_exc()

# Test torch CUDA
print("\nTesting PyTorch CUDA...")
try:
    import torch
    print(f"  CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
except Exception as e:
    print(f"  PyTorch CUDA: FAILED ({e})")

all_ok = all("FAILED" not in str(v) for v in results.values())
print(f"\n{'ALL PACKAGES OK!' if all_ok else 'Some packages had issues'}")
PYEOF

echo ""
echo "=== Complete at $(date) ==="
