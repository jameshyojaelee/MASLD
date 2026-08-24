#!/usr/bin/env bash
#SBATCH --job-name=masld-gse296875-smoke
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PACKAGE_ROOT="${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark"
SOURCE_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration"
DOWNLOAD_ROOT="${SOURCE_ROOT}/candidates/external-gse296875-download-v44-2026-08-18T000000Z"
MAPPING_ROOT="${SOURCE_ROOT}/candidates/external-gse296875-metadata-v45-2026-08-17T230000Z"
LABEL_ROOT="${SOURCE_ROOT}/candidates/external-gse296875-author-labels-v46-2026-08-17T234000Z"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="${PACKAGE_ROOT}/executions/gse296875-rna-atac-smoke-${RUN_STAMP}-${SLURM_JOB_ID}"
WORK_DIR="${RUN_DIR}/work"

mkdir -p "${WORK_DIR}"
printf '%s\n' "${RUN_DIR}" > "${PACKAGE_ROOT}/executions/masld-gse296875-smoke-${SLURM_JOB_ID}.run_dir"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${PACKAGE_ROOT}/src"

"${PYTHON_BIN}" "${PACKAGE_ROOT}/scripts/build_gse296875_rna_atac_smoke.py" select \
  --mapping-metadata "${MAPPING_ROOT}/mapping_metadata.tsv" \
  --author-labels "${LABEL_ROOT}/author_labels.tsv" \
  --output "${WORK_DIR}/selection.tsv"

printf '%s  %s\n' \
  "6af74b7dbdd3839841f90cb7f3e47a3184b4dd3dd02242f6f3bdbd14533cfcab" \
  "${DOWNLOAD_ROOT}/GSE296875_processed.rds" \
  | sha256sum --check --strict

module load R/4.4.3
Rscript "${PACKAGE_ROOT}/scripts/export_gse296875_atac_smoke.R" \
  "${DOWNLOAD_ROOT}/GSE296875_processed.rds" \
  "${WORK_DIR}/selection.tsv" \
  "${WORK_DIR}/atac_export"

"${PYTHON_BIN}" "${PACKAGE_ROOT}/scripts/build_gse296875_rna_atac_smoke.py" assemble \
  --selection "${WORK_DIR}/selection.tsv" \
  --source-lock "${SOURCE_ROOT}/reference/gse296875_source_lock_v44.json" \
  --atac-export "${WORK_DIR}/atac_export" \
  --output "${RUN_DIR}/subset" \
  --authority-root "${RUN_DIR}/proposed_authorities" \
  --generated-contract "${RUN_DIR}/gse296875_rna_atac_smoke_1000_v1.json"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${PACKAGE_ROOT}/src" python -c \
  'import json,sys; from masld_bench.artifacts import verify_frozen_tree; print(json.dumps(verify_frozen_tree(sys.argv[1]), sort_keys=True))' \
  "${RUN_DIR}/subset" \
  > "${RUN_DIR}/control_plane_verification.json"

"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"
sha256sum \
  "${PACKAGE_ROOT}/scripts/build_gse296875_rna_atac_smoke.py" \
  "${PACKAGE_ROOT}/scripts/export_gse296875_atac_smoke.R" \
  "${PACKAGE_ROOT}/slurm/build_gse296875_rna_atac_smoke.sbatch" \
  > "${RUN_DIR}/source_sha256.txt"
