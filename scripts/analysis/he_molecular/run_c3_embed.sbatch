#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=l40s
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
# C3b full: tiles, fat fraction and Phikon-v2 embeddings for all 224 C1-matched slides.
# Usage: sbatch --job-name=<n> --output=<c2>/logs/%j.out run_c3_embed.sbatch <c1_exec_dir> <c2_exec_dir>
set -euo pipefail
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$repo"
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/.local/bin/micromamba run -n he_molecular python -m pip freeze > "$2/pip_freeze_he_molecular.txt"
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/.local/bin/micromamba run -n he_molecular python scripts/analysis/he_molecular/c3_tiles_embed.py \
  --matched "$1/liver_slides_matched.tsv" --seal scripts/analysis/he_molecular/donor_seal.tsv \
  --slides "$repo/data/external/gtex_liver_slides/series" --out "$2/embed"
