#!/usr/bin/env bash
#SBATCH --job-name=scanpy
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/vu_imexcl_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/vu_imexcl_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

python -u Analysis/Spatial/scripts/37_vu_immune_exclusion.py "$@"
