#!/usr/bin/env bash
#SBATCH --job-name=zenodo
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/public_expansion_2026/zenodo_visium_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/public_expansion_2026/zenodo_visium_%j.err

set -euo pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p "${PROJECT_ROOT}/Analysis/Spatial/logs/public_expansion_2026"
bash "${PROJECT_ROOT}/Analysis/Spatial/scripts/public_expansion_2026/download_yakubovsky_visium.sh"
