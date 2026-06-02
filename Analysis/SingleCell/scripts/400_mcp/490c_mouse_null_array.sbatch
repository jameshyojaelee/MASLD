#!/bin/bash
#SBATCH --job-name=mcp_mouse_null_vec
#SBATCH --array=0-9
#SBATCH --partition=gpu
#SBATCH --qos=interactive
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/490c_null_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/490c_null_%A_%a.err
set +u
source "$HOME/.bashrc"
micromamba activate rapids_singlecell
set -u

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

SEED_OFFSET=$(( SLURM_ARRAY_TASK_ID * 100 ))
OUT="${MASLD_PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/mcp/validation/mouse/null_chunks/null_chunk_${SLURM_ARRAY_TASK_ID}.tsv"

echo "[490c.array] task_id=${SLURM_ARRAY_TASK_ID} seed_offset=${SEED_OFFSET} out=${OUT}"
echo "[490c.array] host=$(hostname)"
nvidia-smi -L || true

python 490c_mouse_null_vectorized.py \
    --mode null \
    --k 16 \
    --n-perms 100 \
    --seed-offset "${SEED_OFFSET}" \
    --output "${OUT}" \
    --use-gpu

echo "[490c.array] DONE."
