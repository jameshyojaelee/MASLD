#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --time=72:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/gears_phase2_train_%j.log

# Phase C — production fine-tune on full Saunders Perturb-Multi hep AnnData.
# 20 epochs, 10% holdout (train_gene_set_size=0.9 in prepare_split).

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -u

python Analysis/Perturbation/scripts/d3_synergy/gears_finetune.py \
  --epochs 20 --lr 1e-3 --batch-size 32 \
  --holdout-frac 0.10 \
  --out-dir data/perturbation/finetuned/gears/saunders_hep_v1
