#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/gears_phase2_smoke_%j.log

# Phase B — smoke fine-tune: 1k Saunders cells, 2 epochs.
# Validates the pipeline end-to-end before committing to a long Phase C run.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -u

python Analysis/Perturbation/scripts/d3_synergy/gears_finetune.py \
  --smoke --smoke-cells 1000 \
  --epochs 2 --lr 1e-3 --batch-size 32 \
  --out-dir data/perturbation/finetuned/gears/saunders_hep_smoke
