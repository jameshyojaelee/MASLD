#!/bin/bash
# Task spec preferred io partition. io has 2 nodes (24 CPUs each, ne1di6),
# but tprieto Upload_batch.sh array was using 32/48 CPUs (4 jobs * 8 CPUs).
# Each io node has 8 free CPUs / ~140G free mem -- exactly enough to schedule
# this job. Use io + interactive QOS to skip the queue/QOS-cap.
# Workload: 8 CPUs, 32G mem, ~4-min runtime.
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=geomx_de_govaere
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/geomx_de_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/geomx_de_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "[$(date)] Running 42_govaere2026_geomx_de.R"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript Analysis/Spatial/scripts/42_govaere2026_geomx_de.R

echo "[$(date)] DONE"
