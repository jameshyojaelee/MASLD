#!/usr/bin/env bash
#SBATCH --job-name=mouse_umap_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Integration/mouse_umap_gpu_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Integration/mouse_umap_gpu_%j.err

# GPU-accelerated mouse UMAP integration
# Uses rapids_singlecell for preprocessing + embedding (16-60x speedup)
# Harmony batch correction stays on CPU (fast enough)

set -euo pipefail

echo "=== Mouse UMAP GPU Pipeline ==="
echo "Node: $(hostname) | Date: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

export PYTHONUNBUFFERED=1
export USE_GPU=1

# Add pip-installed NVIDIA libs to LD_LIBRARY_PATH
SITE_PKGS="${CONDA_PREFIX}/lib/python3.12/site-packages"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python3 "${SCRIPT_DIR}/generate_mouse_umap.py"

echo "Done: $(date)"
