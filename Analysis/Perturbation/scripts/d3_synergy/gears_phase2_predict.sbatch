#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/gears_phase2_predict_%j.log

# Phase D — load Saunders-tuned checkpoint, GI_predict on D3 tier1 + tier2 invocab pairs.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -u

CKPT_DIR=${CKPT_DIR:-data/perturbation/finetuned/gears/saunders_hep_v1}
DATASET=${DATASET:-saunders_hep}

# Tier 1
python Analysis/Perturbation/scripts/d3_synergy/gears_predict_v2.py \
  --checkpoint-dir "$CKPT_DIR" \
  --dataset-name "$DATASET" \
  --hits Analysis/Perturbation/data/hits/d3_tier1_pairs_invocab.csv \
  --tier-label tier1_pairs \
  --out-json Analysis/Perturbation/results/d3_synergy/gears_v2_finetune_tier1_pairs__invocab.json

# Tier 2
python Analysis/Perturbation/scripts/d3_synergy/gears_predict_v2.py \
  --checkpoint-dir "$CKPT_DIR" \
  --dataset-name "$DATASET" \
  --hits Analysis/Perturbation/data/hits/d3_tier2_pairs_invocab.csv \
  --tier-label tier2_pairs \
  --out-json Analysis/Perturbation/results/d3_synergy/gears_v2_finetune_tier2_pairs__invocab.json

# Saunders-internal top-200 (Saunders-vocab-compatible pairs; expected non-empty output)
python Analysis/Perturbation/scripts/d3_synergy/gears_predict_v2.py \
  --checkpoint-dir "$CKPT_DIR" \
  --dataset-name "$DATASET" \
  --hits Analysis/Perturbation/data/hits/d3_saunders_internal_pairs_top200.csv \
  --tier-label saunders_internal \
  --out-json Analysis/Perturbation/results/d3_synergy/gears_v2_finetune_saunders_internal__top200.json
