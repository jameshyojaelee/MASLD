#!/usr/bin/env bash
#SBATCH --job-name=masld-gse238219-topology
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
INPUT_ROOT="${ROOT}/executions/gene-perturbation-source-gse238219-21066693_0"
SOURCE="${INPUT_ROOT}/raw/GSE238219_RAW.tar"
SAMPLES=(GSM7660623 GSM7660624 GSM7660625 GSM7660626 GSM7660627)
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
TARGET="${ROOT}/executions/gse238219-topology-${SAMPLE}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite sample topology: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${INPUT_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
print("source frozen-tree verification passed")
PY
PYTHONPATH=src python scripts/audit_gse238219_sample_topology.py \
  --source "${SOURCE}" \
  --sample "${SAMPLE}" \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse238219_sample_topology.py \
  slurm/audit_gse238219_sample_topology.sbatch \
  "${INPUT_ROOT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" "${SAMPLE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse238219_sample_topology",
    "sample": sys.argv[2],
    "expression_values_read": False,
    "cells_as_replicates": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
