#!/usr/bin/env bash
#SBATCH --job-name=cosmx_macsub
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=524288
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/cosmx_macsub_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/cosmx_macsub_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

python -u Analysis/Spatial/scripts/41b_govaere2026_cosmx_mac_subcluster.py "$@"
