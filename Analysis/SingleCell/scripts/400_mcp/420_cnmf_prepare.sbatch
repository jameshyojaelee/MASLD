#!/bin/bash
#SBATCH --job-name=mcp_cnmf_prep
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=300G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_prepare_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_prepare_%j.err

# Arguments (env vars): CNMF_INPUT, CNMF_NAME, CNMF_KLIST, CNMF_NITER, CNMF_NHVG
set -euo pipefail
CNMF_INPUT="${CNMF_INPUT:-atlas_cnmf_global.h5ad}"
CNMF_NAME="${CNMF_NAME:-global}"
CNMF_KLIST="${CNMF_KLIST:-10:12:14:16:18:20:22:25:28:32:36:40}"
CNMF_NITER="${CNMF_NITER:-200}"
CNMF_NHVG="${CNMF_NHVG:-2500}"

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
micromamba activate rapids_singlecell

echo "[420.prepare] host=$(hostname) input=${CNMF_INPUT} name=${CNMF_NAME} k=${CNMF_KLIST}"
python 420_run_cnmf.py prepare \
  --input "${CNMF_INPUT}" \
  --name "${CNMF_NAME}" \
  --k "${CNMF_KLIST}" \
  --n-iter "${CNMF_NITER}" \
  --num-hvg "${CNMF_NHVG}"

echo "[420.prepare] DONE."
