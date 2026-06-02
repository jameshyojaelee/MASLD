#!/bin/bash
#SBATCH --job-name=liana_causal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/liana_causal_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/liana_causal_%j.err

set -euo pipefail

SC_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell"
mkdir -p "${SC_DIR}/logs"

echo "=== LIANA Causal Network Analysis ==="
echo "Start: $(date)"

module purge 2>/dev/null || true
eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell
export PYTHONNOUSERSITE=1

# LD_LIBRARY_PATH fix for pip-installed NVIDIA libs
export LD_LIBRARY_PATH="$(python -c 'import nvidia; import os; print(os.path.dirname(nvidia.__file__))' 2>/dev/null)/cublas/lib:${LD_LIBRARY_PATH:-}" 2>/dev/null || true

python "${SC_DIR}/scripts/fig2_liana_causalnet.py"

echo "=== Done: $(date) ==="
ls -la "${SC_DIR}/results_gpu_v2/fig2_data/liana_"*.csv
