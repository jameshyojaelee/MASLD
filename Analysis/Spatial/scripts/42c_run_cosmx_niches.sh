#!/usr/bin/env bash
# Spatial niche tests on Govaere 2026 CosMx (522,145 cells x 968 genes).
# io partition + interactive QOS — bypasses 7TB nslab cap; 4 CPUs / 64G is enough.
#SBATCH --job-name=cosmx_niche_govaere
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/cosmx_niche_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/cosmx_niche_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

echo "[$(date)] 42c_govaere2026_cosmx_spatial_niches.py"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

python -u Analysis/Spatial/scripts/42c_govaere2026_cosmx_spatial_niches.py "$@"

echo "[$(date)] DONE"
