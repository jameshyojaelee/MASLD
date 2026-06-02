#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --job-name=geomx_load_govaere
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/geomx_load_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/geomx_load_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook --shell bash)"
micromamba activate geomx

echo "[$(date)] Running 40_govaere2026_geomx_load.R"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript Analysis/Spatial/scripts/40_govaere2026_geomx_load.R

echo "[$(date)] DONE"
