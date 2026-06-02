#!/usr/bin/env bash
#SBATCH --job-name=pragmatic_fix
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/pragmatic_fix_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/pragmatic_fix_%j.err
set -euo pipefail

module load CUDA/12.1.1

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

PIP_NVIDIA="${CONDA_PREFIX}/lib/python3.12/site-packages/nvidia"
for pkg in cublas cusparse cusolver curand cufft nvjitlink cusparselt; do
    d="${PIP_NVIDIA}/${pkg}/lib"
    [ -d "$d" ] && export LD_LIBRARY_PATH="${d}:${LD_LIBRARY_PATH:-}"
done

# Pragmatic approach: install squidpy 1.6.x with spatialdata
# Accept dask upgrade from 2024.12.1 → 2025.x
# RAPIDS cugraph/cucim are not used by the spatial pipeline
# The other RAPIDS packages (cuml, rapids_singlecell) will still work

echo "=== Reinstalling squidpy 1.6.x with deps ==="
pip install "squidpy>=1.6,<1.7" spatialdata "dask>=2025.2,<2026" distributed 2>&1

echo ""
echo "=== Verify ALL packages ==="
python << 'PYEOF'
import traceback

results = {}
for name, imp in [("squidpy", "squidpy"), ("cell2location", "cell2location"),
                   ("liana", "liana"), ("hotspot", "hotspot"),
                   ("scanpy", "scanpy"), ("scvi-tools", "scvi"),
                   ("rapids_singlecell", "rapids_singlecell"),
                   ("dask", "dask")]:
    try:
        mod = __import__(imp)
        ver = getattr(mod, "__version__", "OK")
        results[name] = ver
        print(f"  {name}: {ver}")
    except Exception as e:
        results[name] = f"FAILED"
        print(f"  {name}: FAILED ({str(e)[:80]})")
        traceback.print_exc()

# Test squidpy spatial functions
print("\nTesting squidpy.gr functions...")
try:
    from squidpy import gr
    funcs = ["spatial_neighbors", "spatial_autocorr", "nhood_enrichment", "ligrec"]
    for f in funcs:
        assert hasattr(gr, f), f"{f} missing"
    print("  All squidpy.gr spatial analysis functions: OK")
except Exception as e:
    print(f"  squidpy.gr: FAILED ({e})")
    traceback.print_exc()

# Test rapids_singlecell basic functions
print("\nTesting rapids_singlecell...")
try:
    import rapids_singlecell as rsc
    from rapids_singlecell import pp, tl
    print("  rapids_singlecell pp/tl: OK")
except Exception as e:
    print(f"  rapids_singlecell: FAILED ({e})")

# Test cell2location model import
print("\nTesting cell2location model import...")
try:
    from cell2location.models import RegressionModel, Cell2location
    print("  cell2location models: OK")
except Exception as e:
    print(f"  cell2location models: FAILED ({e})")

all_ok = all("FAILED" not in str(v) for v in results.values())
print(f"\n{'ALL PACKAGES OK!' if all_ok else 'Some packages had issues'}")
PYEOF

echo ""
echo "=== Complete at $(date) ==="
