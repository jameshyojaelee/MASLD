#!/bin/bash
#SBATCH --job-name=spatial_zonation_rerun
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/spatial_zonation_rerun_%j.out
#SBATCH --error=logs/spatial_zonation_rerun_%j.err

# Re-runs zonation pipeline after GLS2 was moved from pericentral_markers
# to periportal_markers in Analysis/Spatial/config/spatial_params.yaml.
# Runs: 04a (define zonation) -> 05b (ligand-receptor) -> 09a (zonation-aware DE)
# NOTE: Scripts 04a/05b/09a do not support --force; they overwrite outputs unconditionally.

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial

echo "=== Step 1: Define liver zonation axis (04a) ==="
echo "Started: $(date)"
python scripts/04a_define_zonation.py 2>&1
echo "Completed 04a: $(date)"

echo "=== Step 2: Ligand-receptor interaction analysis (05b) ==="
echo "Started: $(date)"
python scripts/05b_ligand_receptor.py 2>&1
echo "Completed 05b: $(date)"

echo "=== Step 3: Zonation-aware differential expression (09a) ==="
echo "Started: $(date)"
python scripts/09a_zonation_aware_de.py 2>&1
echo "Completed 09a: $(date)"

echo "=== Spatial zonation rerun complete: $(date) ==="
