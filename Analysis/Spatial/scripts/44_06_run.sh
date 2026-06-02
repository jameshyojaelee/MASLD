#!/bin/bash
#SBATCH --job-name=spatial_govaere_integ
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/spatial_govaere_integ_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/spatial_govaere_integ_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== 44_spatial_consensus ==="
python Analysis/Spatial/scripts/44_spatial_consensus.py

echo "=== 06_integration ==="
python Analysis/Spatial/scripts/06_integration.py

echo "=== DONE ==="
