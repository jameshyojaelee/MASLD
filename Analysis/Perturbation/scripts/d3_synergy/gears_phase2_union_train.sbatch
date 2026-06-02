#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --time=90:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/gears_phase2_union_train_%j.log

# Phase C — production 3-way union-vocab fine-tune
# (Replogle K562 ∪ Norman K562 ∪ Saunders mouse hep; full data; 50 epochs).
# Output: data/perturbation/finetuned/gears/saunders_hep_v1/best.ckpt

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -u

python Analysis/Perturbation/scripts/d3_synergy/gears_union_finetune.py \
  --epochs 50 --lr 1e-3 --batch-size 32 \
  --holdout-frac 0.10 \
  --out-dir data/perturbation/finetuned/gears/saunders_hep_v1
