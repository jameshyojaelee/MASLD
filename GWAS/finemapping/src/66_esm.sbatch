#!/bin/bash
#SBATCH --job-name=esm
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/esm_%j.log
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
# fair-esm staged (--no-deps) alongside rapids_singlecell's torch; offline weights
export PYTHONPATH="$ROOT/.envs/esm_pkg:${PYTHONPATH:-}"
export TORCH_HOME="$ROOT/data/external/esm_weights/torch_home"

echo "[esm] host=$(hostname) gpu=${CUDA_VISIBLE_DEVICES:-none}"
micromamba run -n rapids_singlecell python "$ROOT/GWAS/finemapping/src/66_coding_hardening.py" --mode esm
echo "[esm] done"
