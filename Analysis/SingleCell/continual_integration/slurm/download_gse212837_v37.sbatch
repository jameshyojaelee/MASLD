#!/usr/bin/env bash
#SBATCH --job-name=masld-gse212837
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/gse212837-v37-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/gse212837-v37-%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PIPELINE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
OUTPUT="${PIPELINE_ROOT}/candidates/external-gse212837-download-v37-2026-08-15T160000Z"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212837/suppl/GSE212837_all_data_labeled_processed.h5ad.gz"

test ! -e "${OUTPUT}"
mkdir "${OUTPUT}"
export PYTHONNOUSERSITE=1 PYTHONPATH="${PIPELINE_ROOT}" PYTHONUNBUFFERED=1

curl --fail --location --retry 10 --retry-all-errors \
  --output "${OUTPUT}/GSE212837_all_data_labeled_processed.h5ad.gz" "${URL}"
gzip --test "${OUTPUT}/GSE212837_all_data_labeled_processed.h5ad.gz"
sha256sum "${OUTPUT}/GSE212837_all_data_labeled_processed.h5ad.gz" \
  > "${OUTPUT}/download_sha256.txt"
gzip --decompress --stdout "${OUTPUT}/GSE212837_all_data_labeled_processed.h5ad.gz" \
  > "${OUTPUT}/GSE212837_all_data_labeled_processed.h5ad"
sha256sum "${OUTPUT}/GSE212837_all_data_labeled_processed.h5ad" \
  > "${OUTPUT}/h5ad_sha256.txt"
"${PYTHON_BIN}" "${PIPELINE_ROOT}/scripts/inspect_external_gse212837_v37.py" \
  --h5ad "${OUTPUT}/GSE212837_all_data_labeled_processed.h5ad" \
  > "${OUTPUT}/h5ad_inspection.json"
"${PYTHON_BIN}" -m pip freeze > "${OUTPUT}/pip_freeze.txt"
sha256sum "$0" "${PIPELINE_ROOT}/scripts/inspect_external_gse212837_v37.py" \
  > "${OUTPUT}/source_sha256.txt"
touch "${OUTPUT}/DOWNLOAD_AND_INSPECTION_COMPLETE"
