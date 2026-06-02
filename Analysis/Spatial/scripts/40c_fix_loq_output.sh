#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=geomx_loq_fix
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/geomx_loq_fix_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/geomx_loq_fix_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook --shell bash)"
micromamba activate geomx

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript Analysis/Spatial/scripts/40_govaere2026_geomx_load.R
echo "[$(date)] DONE"
