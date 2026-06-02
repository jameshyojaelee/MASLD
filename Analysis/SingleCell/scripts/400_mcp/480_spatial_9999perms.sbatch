#!/bin/bash
#SBATCH --job-name=mcp_spatial_9999
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/480_spatial_9999_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/480_spatial_9999_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"
source "$HOME/.bashrc"
micromamba activate spatial || micromamba activate rapids_singlecell

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "[480.9999] host=$(hostname) k=16 permutations=9999"
python 480_spatial_validation.py --name global --k 16 --top-n 100 --permutations 9999

# Echo expects this filename per mission spec
SRC="${MASLD_PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/mcp/validation/spatial/spatial_moran_global_k16_9999perms.tsv"
DEST="${MASLD_PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/mcp/validation/spatial/spatial_moran_global_k16_9999perms.tsv"
if [ -f "$SRC" ]; then echo "[480.9999] output at $SRC"; fi
echo "[480.9999] DONE."
