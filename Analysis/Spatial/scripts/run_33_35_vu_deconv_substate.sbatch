#!/usr/bin/env bash
#SBATCH --job-name=vu_m3_deconv
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=interactive
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/vu_m3_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/vu_m3_%j.log
set -euo pipefail

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs

module load CUDA/12.1.1 || true

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial

echo "=== Team M3: Vu cell-type composition decoupling ==="
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'no-gpu')"

# Step 1: Cell2location deconvolution on Vu spots
echo ""
echo "--- 33 cell2location deconvolution ---"
python scripts/33_vu_cell2location_deconvolution.py --max-epochs 15000 \
    || python scripts/33_vu_cell2location_deconvolution.py --use-fallback

# Step 2: Within-hepatocyte bimodality (CPU)
echo ""
echo "--- 34 within-hepatocyte bimodality ---"
python scripts/34_vu_substate_within_hepatocytes.py --n-mc 999

# Step 3: Matched-null Moran's I (CPU)
echo ""
echo "--- 35 matched-null Moran's I ---"
python scripts/35_vu_moran_matched_null.py --n-null 100

# Step 4: Aggregate verdict
echo ""
echo "--- 36 summary verdict ---"
python scripts/36_vu_substate_summary.py

echo ""
echo "Finished: $(date)"
