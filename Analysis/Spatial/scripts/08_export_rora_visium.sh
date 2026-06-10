#!/bin/bash
#SBATCH --job-name=roravis
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/roravis_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts/logs/roravis_%j.err

set -eo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$BASE/Analysis/Spatial/scripts/logs"
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial
python "$BASE/Analysis/Spatial/scripts/08_export_rora_visium.py"
