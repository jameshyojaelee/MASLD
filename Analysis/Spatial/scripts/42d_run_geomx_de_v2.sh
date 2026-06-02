#!/bin/bash
# v2 GeoMx DE: limma+dupCor, permFDR, 1-sided canonical markers.
# io partition + interactive QOS, 8 CPUs, 32G mem.
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --job-name=geomx_de_v2
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/geomx_de_v2_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/geomx_de_v2_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "[$(date)] Running 42d_govaere2026_geomx_de_v2.R"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript Analysis/Spatial/scripts/42d_govaere2026_geomx_de_v2.R

echo "[$(date)] DONE"
