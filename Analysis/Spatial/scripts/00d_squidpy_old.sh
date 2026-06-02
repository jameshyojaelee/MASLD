#!/usr/bin/env bash
#SBATCH --job-name=squidpy_old
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/squidpy_old_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/squidpy_old_%j.err
set -euo pipefail

module load CUDA/12.1.1

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

PIP_NVIDIA="${CONDA_PREFIX}/lib/python3.12/site-packages/nvidia"
for pkg in cublas cusparse cusolver curand cufft nvjitlink cusparselt; do
    d="${PIP_NVIDIA}/${pkg}/lib"
    [ -d "$d" ] && export LD_LIBRARY_PATH="${d}:${LD_LIBRARY_PATH:-}"
done

# Try progressively older squidpy versions until one works
for VER in "1.5.0" "1.4.1" "1.4.0" "1.3.1" "1.3.0" "1.2.4" "1.2.3"; do
    echo "=== Trying squidpy ${VER} ==="
    pip install "squidpy==${VER}" --force-reinstall --no-deps 2>&1 || continue

    python -c "
import squidpy as sq
print(f'squidpy {sq.__version__} imports OK')
from squidpy import gr
assert hasattr(gr, 'spatial_neighbors')
assert hasattr(gr, 'spatial_autocorr')
assert hasattr(gr, 'nhood_enrichment')
assert hasattr(gr, 'ligrec')
print('All squidpy.gr functions available')
" 2>&1 && echo "SUCCESS: squidpy ${VER} works!" && break || echo "FAILED: squidpy ${VER}"
    echo ""
done

echo ""
echo "=== Final verification ==="
python << 'PYEOF'
for name, imp in [("squidpy", "squidpy"), ("cell2location", "cell2location"),
                   ("liana", "liana"), ("hotspot", "hotspot"),
                   ("rapids_singlecell", "rapids_singlecell"),
                   ("dask", "dask"), ("scanpy", "scanpy")]:
    try:
        mod = __import__(imp)
        print(f"  {name}: {getattr(mod, '__version__', 'OK')}")
    except Exception as e:
        print(f"  {name}: FAILED ({str(e)[:60]})")
PYEOF

echo "=== Complete at $(date) ==="
