#!/usr/bin/env bash
#SBATCH --job-name=diag_dask
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/diag_dask_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/diag_dask_%j.err
set -euo pipefail

module load CUDA/12.1.1

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

PIP_NVIDIA="${CONDA_PREFIX}/lib/python3.12/site-packages/nvidia"
for pkg in cublas cusparse cusolver curand cufft nvjitlink cusparselt; do
    d="${PIP_NVIDIA}/${pkg}/lib"
    [ -d "$d" ] && export LD_LIBRARY_PATH="${d}:${LD_LIBRARY_PATH:-}"
done

echo "=== Package versions ==="
python -c "import dask; print(f'dask: {dask.__version__}')"
python -c "
try:
    import dask_expr
    print(f'dask_expr: INSTALLED')
except ImportError:
    print('dask_expr: NOT INSTALLED')
"

echo ""
echo "=== Full traceback: import dask.dataframe ==="
python -c "
import traceback
try:
    import dask.dataframe as dd
    print('dask.dataframe: OK')
    print(f'  type: {type(dd.DataFrame)}')
except Exception as e:
    print(f'dask.dataframe: FAILED')
    traceback.print_exc()
"

echo ""
echo "=== Full traceback: import rapids_singlecell ==="
python -c "
import traceback
try:
    import rapids_singlecell
    print(f'rapids_singlecell: {rapids_singlecell.__version__}')
except Exception as e:
    print(f'rapids_singlecell: FAILED')
    traceback.print_exc()
"

echo ""
echo "=== Full traceback: import squidpy ==="
python -c "
import traceback
try:
    import squidpy
    print(f'squidpy: {squidpy.__version__}')
except Exception as e:
    print(f'squidpy: FAILED')
    traceback.print_exc()
"

echo ""
echo "=== Check dask._dask_expr_enabled ==="
python -c "
import dask
print(f'dask version: {dask.__version__}')
print(dir(dask))
if hasattr(dask, '_dask_expr_enabled'):
    print(f'_dask_expr_enabled: {dask._dask_expr_enabled()}')
"

echo ""
echo "=== Check dask legacy/expr switch ==="
python -c "
import dask
# Check if there's a way to use legacy dask
try:
    import dask.dataframe._compat
    print('_compat module exists')
except:
    pass
try:
    dask.config.set({'dataframe.query-planning': False})
    print('Set query-planning=False')
    import dask.dataframe as dd
    print(f'dask.dataframe with query-planning=False: OK')
except Exception as e:
    print(f'Failed: {e}')
"

echo ""
echo "=== Complete at $(date) ==="
