#!/bin/bash
#SBATCH --job-name=polars
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.log
#
# End-to-end Cas13 guide build:
#   1. build the parquet index from the upstream genome-wide guide CSV (one-time)
#   2. build the full library + control guide tables
#
# The eager sorted index build holds the slim ~12.6M-row table in memory, hence
# 120 G. Run with:  sbatch build_guides.sbatch [RELEASE]   (default vM38)
# NB: no `-u` — micromamba/conda activation scripts reference unbound vars.
set -eo pipefail

RELEASE="${1:-vM38}"
# SLURM runs a spool copy of this script, so BASH_SOURCE is unreliable; pin the
# scripts dir by absolute path (MASLD_PROJECT_ROOT override honoured).
HERE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}/Cas13_Library_Design/scripts/guides"
cd "$HERE"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== [1/2] build_guide_index ($RELEASE) ==="
python build_guide_index.py --release "$RELEASE"

echo "=== [2/2] build_library_guides ($RELEASE) ==="
python build_library_guides.py --release "$RELEASE"

echo "=== ALL DONE ==="
