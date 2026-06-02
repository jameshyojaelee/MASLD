#!/bin/bash
#SBATCH --job-name=spatial_rebuild
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/spatial_rebuild_%j.out
#SBATCH --error=logs/spatial_rebuild_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== Step 1: Rebuild AnnData (HVG fix) ==="
echo "Started: $(date)"
python scripts/02_build_anndata.py --force 2>&1
echo "Completed 02: $(date)"

echo "=== Step 2: Prepare cell2location reference (stratified subsampling fix) ==="
echo "Started: $(date)"
python scripts/03a_prepare_c2l_reference.py --force 2>&1
echo "Completed 03a: $(date)"

echo "Spatial rebuild complete: $(date)"
