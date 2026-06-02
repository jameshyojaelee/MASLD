#!/bin/bash
#SBATCH --job-name=test_imports
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/test_import_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/test_import_%j.err

set -euo pipefail

CONDA_ENV="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate "${CONDA_ENV}"
export LD_LIBRARY_PATH="${CONDA_ENV}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

echo "Testing imports on $(hostname)..."
echo "CPU: $(lscpu | grep 'Model name')"

python3 -c "
import sys; print(f'Python: {sys.version}')
print('Importing numpy...'); import numpy; print(f'  numpy {numpy.__version__}')
print('Importing scanpy...'); import scanpy; print(f'  scanpy {scanpy.__version__}')
print('Importing scvi...'); import scvi; print(f'  scvi {scvi.__version__}')
print('Importing celltypist...'); import celltypist; print(f'  celltypist {celltypist.__version__}')
print('All imports OK')
"
