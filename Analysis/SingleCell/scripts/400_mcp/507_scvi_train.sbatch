#!/bin/bash
#SBATCH --job-name=mcp_scvi_train
#SBATCH --partition=gpu
#SBATCH --qos=nslab
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/507_scvi_train_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/507_scvi_train_%j.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

# Add pip-installed NVIDIA libs to LD_LIBRARY_PATH (required for torch/scvi on GPU)
SITE_PKGS="${CONDA_PREFIX}/lib/python3.12/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "[507] host=$(hostname)"
nvidia-smi || true

python 507_scvi_subsample_cnmf.py --stage train --n-target 200000 --n-latent 30 --max-epochs 400 --seed 42

echo "[507] DONE — denoised h5ad ready for cNMF chain."
