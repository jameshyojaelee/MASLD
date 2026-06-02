#!/bin/bash
#SBATCH --job-name=mcp_cnmf_cons
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_consensus_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_consensus_%j.err

# Env: CNMF_NAME (required), CNMF_KLIST, CNMF_DTHRESH
set -euo pipefail
CNMF_NAME="${CNMF_NAME:-global}"
CNMF_KLIST="${CNMF_KLIST:-10:12:14:16:18:20:22:25:28:32:36:40}"
CNMF_DTHRESH="${CNMF_DTHRESH:-0.03}"

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
micromamba activate rapids_singlecell

python 420_run_cnmf.py consensus \
  --name "${CNMF_NAME}" \
  --k "${CNMF_KLIST}" \
  --density-threshold "${CNMF_DTHRESH}"

echo "[420.consensus] DONE."
