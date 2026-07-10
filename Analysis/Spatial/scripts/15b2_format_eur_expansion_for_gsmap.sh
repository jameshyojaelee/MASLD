#!/usr/bin/env bash
#SBATCH --job-name=gsmap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/gsmap_fmt_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/gsmap_fmt_%j.err
set -eo pipefail
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
PR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PR"
python "$PR/Analysis/Spatial/scripts/15b2_format_eur_expansion_for_gsmap.py"
