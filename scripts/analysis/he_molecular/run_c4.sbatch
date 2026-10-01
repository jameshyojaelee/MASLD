#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
# C4: Model C development analysis (ridge, shuffle null, BH). Reads development files only.
# Usage: sbatch --job-name=<n> --output=<c2>/logs/%j.out run_c4.sbatch <c2_exec_dir>
set -euo pipefail
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/.local/bin/micromamba run -n he_molecular python \
  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/he_molecular/c4_ridge.py \
  --targets "$1/targets/targets_development.tsv" --ceilings "$1/targets/program_ceilings.tsv" \
  --features "$1/embed/slide_features.tsv.gz" --out "$1/${2:-c4}"
