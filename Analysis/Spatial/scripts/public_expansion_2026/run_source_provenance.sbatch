#!/usr/bin/env bash
#SBATCH --job-name=provenance
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/public_expansion_2026/provenance_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/public_expansion_2026/provenance_%j.err

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"

SPATIAL_PYTHON="${SPATIAL_PYTHON:-/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python}"
test -x "${SPATIAL_PYTHON}" || {
  echo "Missing spatial-environment Python: ${SPATIAL_PYTHON}" >&2
  exit 2
}
"${SPATIAL_PYTHON}" \
  "${PROJECT_ROOT}/Analysis/Spatial/scripts/public_expansion_2026/finalize_source_provenance.py"
