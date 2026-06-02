#!/bin/bash
#SBATCH --job-name=fig3_advanced
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=8:00:00
#SBATCH --output=logs/fig3_advanced_%j.out
#SBATCH --error=logs/fig3_advanced_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell
mkdir -p logs

# Fix LD_LIBRARY_PATH for pip-installed NVIDIA libs
export LD_LIBRARY_PATH="/gpfs/commons/home/jameslee/mambaforge/envs/rapids_singlecell/lib/python3.12/site-packages/nvidia/cusparselt/lib:${LD_LIBRARY_PATH:-}"

echo "=== Fig3 Advanced Analyses ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
nvidia-smi || true

micromamba run -n rapids_singlecell python scripts/fig2_advanced_analyses.py

echo "=== Done: $(date) ==="
ls -la results_gpu_v2/fig2_data/
