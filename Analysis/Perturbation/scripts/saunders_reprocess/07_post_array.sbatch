#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Perturbation/scripts/saunders_reprocess/logs/07_post_array_%j.log

# After 32-SRR demux array completes, aggregate per-cell counts and validate.

set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
WORK=$PROJ/Analysis/Perturbation/scripts/saunders_reprocess
DATA=$PROJ/data/perturbation/datasets/saunders2025

echo "[$(date +%T)] === Aggregating per-cell sgRNA assignments ==="
mkdir -p "$DATA/processed"
python3 "$WORK/04_aggregate_cells.py" \
    --in-dir "$DATA/processed/per_srr_counts" \
    --out "$DATA/processed/per_cell_sgrna.tsv.gz" \
    --min-umis 3 \
    --min-dominant-frac 0.5

echo "[$(date +%T)] === SRR-to-sublib mapping ==="
python3 "$WORK/05_map_srr_to_sublib.py" \
    --counts-dir "$DATA/processed/per_srr_counts" \
    --out "$DATA/processed/srr_sublib_map.json"

echo "[$(date +%T)] === Validating final assignments ==="
python3 "$WORK/06_validate_final.py" \
    --per-cell "$DATA/processed/per_cell_sgrna.tsv.gz" \
    --report-out "$DATA/processed/validation_report.json"

echo "[$(date +%T)] === Final outputs ==="
ls -la "$DATA/processed/"

echo "[$(date +%T)] DONE"
