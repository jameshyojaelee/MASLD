#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/gears_phase2_union_smoke_%j.log

# Phase B — union-vocab smoke fine-tune (Norman ∪ Saunders subset, 2 epochs).
# Validates end-to-end pipeline + union vocab behavior before Phase C.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -u

python Analysis/Perturbation/scripts/d3_synergy/gears_union_finetune.py \
  --smoke --smoke-cells 1500 \
  --epochs 2 --lr 1e-3 --batch-size 32 \
  --out-dir data/perturbation/finetuned/gears/saunders_hep_smoke
