#!/usr/bin/env bash
#SBATCH --job-name=vu_m3_fast
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/vu_m3_fast_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/vu_m3_fast_%j.log
set -euo pipefail

# FAST PATH: signature-based per-spot cell-type scoring instead of cell2location.
# Used as the primary deliverable for Team M3 because cell2location requires
# 12-24 h on GPU. If/when GPU c2l completes, vu_deconvolved_substate.csv can be
# rerun with --force on the GPU job and downstream scripts re-run.

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial

echo "=== Team M3 (FAST signature path): Vu cell-type composition decoupling ==="
echo "Started: $(date)"
echo "Host: $(hostname)"

echo ""
echo "--- 33 signature-based per-spot scoring (fallback) ---"
python scripts/33_vu_cell2location_deconvolution.py --use-fallback --force

echo ""
echo "--- 34 within-hepatocyte bimodality ---"
python scripts/34_vu_substate_within_hepatocytes.py --n-mc 999 --force

echo ""
echo "--- 35 matched-null Moran's I ---"
python scripts/35_vu_moran_matched_null.py --n-null 100 --force

echo ""
echo "--- 36 summary verdict ---"
python scripts/36_vu_substate_summary.py

echo ""
echo "Finished: $(date)"
