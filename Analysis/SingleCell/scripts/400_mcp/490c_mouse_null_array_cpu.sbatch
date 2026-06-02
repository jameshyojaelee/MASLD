#!/bin/bash
# Fallback: pure-CPU array, used if GPU partition is full.
# Vectorized matmul on CPU is still ~50× faster than the 490b score_genes loop.
#SBATCH --job-name=mcp_mouse_null_vec_cpu
#SBATCH --array=0-9
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=210G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/490c_null_cpu_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/490c_null_cpu_%A_%a.err
set +u
source "$HOME/.bashrc"
micromamba activate rapids_singlecell
set -u

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

SEED_OFFSET=$(( SLURM_ARRAY_TASK_ID * 100 ))
OUT="${MASLD_PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/mcp/validation/mouse/null_chunks/null_chunk_cpu_${SLURM_ARRAY_TASK_ID}.tsv"

echo "[490c.array.cpu] task_id=${SLURM_ARRAY_TASK_ID} seed_offset=${SEED_OFFSET} out=${OUT}"
echo "[490c.array.cpu] host=$(hostname)"

python 490c_mouse_null_vectorized.py \
    --mode null \
    --k 16 \
    --n-perms 100 \
    --seed-offset "${SEED_OFFSET}" \
    --output "${OUT}"

echo "[490c.array.cpu] DONE."
