#!/bin/bash
#SBATCH --job-name=mcp_cnmf_fact
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=60G
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --array=0-31
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_factorize_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_factorize_%A_%a.err

# Env: CNMF_NAME (required), TOTAL_WORKERS (default 32; must match --array size)
set -euo pipefail
CNMF_NAME="${CNMF_NAME:-global}"
TOTAL_WORKERS="${TOTAL_WORKERS:-32}"

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
micromamba activate rapids_singlecell
# Prevent accidental thread oversubscription (cNMF uses sklearn which honors these)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

echo "[420.factorize] name=${CNMF_NAME} worker=${SLURM_ARRAY_TASK_ID}/${TOTAL_WORKERS} host=$(hostname)"
python 420_run_cnmf.py factorize \
  --name "${CNMF_NAME}" \
  --worker-index "${SLURM_ARRAY_TASK_ID}" \
  --total-workers "${TOTAL_WORKERS}"

echo "[420.factorize] worker ${SLURM_ARRAY_TASK_ID} DONE."
