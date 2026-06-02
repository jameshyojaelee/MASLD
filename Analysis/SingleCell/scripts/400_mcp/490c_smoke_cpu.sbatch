#!/bin/bash
#SBATCH --job-name=mcp_mouse_smoke_cpu
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=160G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/490c_smoke_cpu_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/490c_smoke_cpu_%j.err
set +u
source "$HOME/.bashrc"
micromamba activate rapids_singlecell
set -u

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

OUT_OBS="${MASLD_PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/mcp/validation/mouse/mouse_observed_directmean_k16.tsv"
OUT_SMOKE="${MASLD_PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/mcp/validation/mouse/null_chunks/null_chunk_smoke_cpu.tsv"

echo "[smoke] computing observed direct-mean (no GPU)"
python 490c_mouse_null_vectorized.py --mode observed --k 16

echo "[smoke] running 3-perm null on CPU"
python 490c_mouse_null_vectorized.py \
    --mode null \
    --k 16 \
    --n-perms 3 \
    --seed-offset 9000 \
    --output "${OUT_SMOKE}"

echo "[smoke] DONE. observed→${OUT_OBS}; null_chunk→${OUT_SMOKE}"
