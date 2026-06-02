#!/bin/bash
#SBATCH --job-name=cnmf_seed2024_fact
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=60G
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --array=0-31
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_seed2024_fact_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/420_seed2024_fact_%A_%a.err
set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"
source "$HOME/.bashrc"
micromamba activate rapids_singlecell
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

python 420_run_cnmf.py factorize \
  --name global_seed2024 \
  --worker-index "${SLURM_ARRAY_TASK_ID}" \
  --total-workers 32

echo "[420.factorize seed2024] worker ${SLURM_ARRAY_TASK_ID} DONE."
