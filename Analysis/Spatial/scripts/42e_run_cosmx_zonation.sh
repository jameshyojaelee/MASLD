#!/usr/bin/env bash
# Hepatocyte zonation scoring on Govaere 2026 CosMx (297K hepatocytes / 522K cells).
# io partition + interactive QOS — bypasses 7TB nslab cap; 4 CPUs, 64G is
# sufficient for in-memory scanpy + KDTree on this dataset.
#SBATCH --job-name=cosmx_zonation
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/cosmx_zonation_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/cosmx_zonation_%j.err

set -eo pipefail
set +u
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

echo "[$(date)] 42e_govaere2026_cosmx_zonation.py"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

python -u Analysis/Spatial/scripts/42e_govaere2026_cosmx_zonation.py "$@"

echo "[$(date)] DONE"
