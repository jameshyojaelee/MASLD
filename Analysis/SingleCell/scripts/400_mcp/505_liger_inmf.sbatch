#!/bin/bash
#SBATCH --job-name=alpha5_pyliger
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/505_pyliger_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/505_pyliger_%j.err
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# Re-stage per-dataset HVG mtx (idempotent; fast; writes liger_staging/)
if [ ! -f /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/mcp/liger_staging/manifest.json ]; then
  python -u 505_liger_prep.py
fi

# pyliger (+ louvain/python-igraph/cmake) pre-installed via conda-forge on the
# login node to avoid source builds. Confirm with the umap.plot stub (same as
# in the main python script) so the import chain doesn't hit cudf.
python -c "
import sys, types
_up = types.ModuleType('umap.plot')
_up.points = lambda *a, **k: None
_up.connectivity = lambda *a, **k: None
_up.diagnostic = lambda *a, **k: None
sys.modules['umap.plot'] = _up
import pyliger
print('[505.sbatch] pyliger import OK')
"

: "${LIGER_K_LIST:=10,16}"
: "${LIGER_N_ITERS:=30}"
export LIGER_K_LIST LIGER_N_ITERS

python -u 505_pyliger_inmf.py
