#!/usr/bin/env bash
#SBATCH --job-name=rm_dask_expr
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/rm_dask_expr_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/rm_dask_expr_%j.err
set -euo pipefail

module load CUDA/12.1.1

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

PIP_NVIDIA="${CONDA_PREFIX}/lib/python3.12/site-packages/nvidia"
for pkg in cublas cusparse cusolver curand cufft nvjitlink cusparselt; do
    d="${PIP_NVIDIA}/${pkg}/lib"
    [ -d "$d" ] && export LD_LIBRARY_PATH="${d}:${LD_LIBRARY_PATH:-}"
done

echo "=== Remove dask-expr (blocks non-2024.12.1 dask) ==="
pip uninstall -y dask-expr 2>&1 || true

echo ""
echo "=== Also update dask-expr-compatible packages ==="
pip install "dask>=2025.2,<2025.6" "distributed>=2025.2,<2025.6" 2>&1

echo ""
echo "=== Verify ==="
python << 'PYEOF'
for name, imp in [("squidpy", "squidpy"), ("cell2location", "cell2location"),
                   ("liana", "liana"), ("hotspot", "hotspot"),
                   ("rapids_singlecell", "rapids_singlecell"),
                   ("dask", "dask"), ("scanpy", "scanpy"), ("scvi", "scvi")]:
    try:
        mod = __import__(imp)
        ver = getattr(mod, "__version__", "OK")
        print(f"  {name}: {ver}")
    except Exception as e:
        print(f"  {name}: FAILED ({str(e)[:80]})")

print("\nTesting squidpy.gr functions...")
try:
    from squidpy import gr
    for f in ["spatial_neighbors", "spatial_autocorr", "nhood_enrichment", "ligrec"]:
        assert hasattr(gr, f), f"{f} missing"
    print("  All squidpy.gr spatial functions: OK")
except Exception as e:
    print(f"  squidpy.gr: FAILED ({str(e)[:80]})")

print("\nTesting rapids_singlecell...")
try:
    import rapids_singlecell as rsc
    from rapids_singlecell import pp, tl
    print("  rapids_singlecell pp/tl: OK")
except Exception as e:
    print(f"  rsc: FAILED ({str(e)[:80]})")

print("\nTesting cell2location models...")
try:
    from cell2location.models import RegressionModel, Cell2location
    print("  cell2location models: OK")
except Exception as e:
    print(f"  c2l: FAILED ({str(e)[:80]})")
PYEOF

echo ""
echo "=== Complete at $(date) ==="
