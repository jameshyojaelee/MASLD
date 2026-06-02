#!/bin/bash
#SBATCH --job-name=mcp_cnmf_pct
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16
#SBATCH --time=90:00:00
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/421_pct_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/421_pct_%A_%a.err

# Per-cell-type cNMF: prepare + internal factorize (single-node, 16 threads) + consensus.
# Array index -> cell-type.
set -euo pipefail

CELLTYPES=(hepatocytes endothelial_cells fibroblasts macrophages cholangiocytes)
CT="${CELLTYPES[$SLURM_ARRAY_TASK_ID]}"
KLIST="${CNMF_KLIST:-5,6,7,8,10,12,15,18,20}"
NITER="${CNMF_NITER:-150}"
NHVG="${CNMF_NHVG:-2000}"

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
micromamba activate rapids_singlecell

NAME="pct_${CT}"
INPUT="atlas_cnmf_${CT}.h5ad"
echo "[421] host=$(hostname) CT=${CT} name=${NAME} k=${KLIST}"

# Prepare
python 420_run_cnmf.py prepare \
  --input "${INPUT}" --name "${NAME}" \
  --k "${KLIST}" --n-iter "${NITER}" --num-hvg "${NHVG}"

# Factorize all replicates on a single 16c node (small enough for per-CT)
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
python 420_run_cnmf.py factorize --name "${NAME}" --worker-index 0 --total-workers 1

# Consensus
python 420_run_cnmf.py consensus --name "${NAME}" --k "${KLIST}" --density-threshold 0.03

echo "[421] ${CT} DONE."
