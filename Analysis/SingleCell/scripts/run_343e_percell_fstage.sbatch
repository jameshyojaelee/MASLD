#!/bin/bash
#SBATCH --job-name=343e_percell_fstage
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343e_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs/343e_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$PROJECT_ROOT"
SCRIPT="$PROJECT_ROOT/Analysis/SingleCell/scripts/343e_percell_fstage_breakdown.py"

echo "[host] $(hostname)"
echo "[date] $(date)"
mkdir -p "$PROJECT_ROOT/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/logs"

export MAMBA_EXE='/gpfs/commons/home/jameslee/.local/bin/micromamba'
export MAMBA_ROOT_PREFIX='/gpfs/commons/home/jameslee/micromamba'
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate rapids_singlecell

python --version
python -c "import sklearn, anndata, numpy, pandas, matplotlib; print('sklearn', sklearn.__version__)"

python "$SCRIPT"

echo "[done] $(date)"
