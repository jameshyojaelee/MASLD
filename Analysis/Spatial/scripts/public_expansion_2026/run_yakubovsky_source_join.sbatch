#!/usr/bin/env bash
#SBATCH --job-name=h5py
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/public_expansion_2026/h5py_join_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/public_expansion_2026/h5py_join_%j.err

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
DATA_DIR="${PROJECT_ROOT}/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026"
SOURCE_DIR="${DATA_DIR}/Human-liver"
OUTPUT_DIR="${YAK_SOURCE_OUTPUT_DIR:-${PROJECT_ROOT}/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/acquisition/yakubovsky2026}"

mkdir -p "${PROJECT_ROOT}/Analysis/Spatial/logs/public_expansion_2026"

/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python \
  "${PROJECT_ROOT}/Analysis/Spatial/scripts/public_expansion_2026/audit_yakubovsky_source_join.py" \
  --data-dir "${DATA_DIR}" \
  --source-dir "${SOURCE_DIR}" \
  --output-dir "${OUTPUT_DIR}"
