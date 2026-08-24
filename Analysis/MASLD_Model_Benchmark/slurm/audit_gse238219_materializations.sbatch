#!/usr/bin/env bash
#SBATCH --job-name=masld-gse238219-cross-sample
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${UPSTREAM_ARRAY_JOB_ID:?UPSTREAM_ARRAY_JOB_ID is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/gse238219-cross-sample-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite cross-sample gate: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python scripts/audit_gse238219_materializations.py \
  --execution-root "${ROOT}/executions" \
  --array-job-id "${UPSTREAM_ARRAY_JOB_ID}" \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse238219_materializations.py \
  slurm/audit_gse238219_materializations.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" "${UPSTREAM_ARRAY_JOB_ID}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse238219_cross_sample_gate",
    "upstream_array_job_id": sys.argv[2],
    "biological_independence_resolved": False,
    "source_perturbation_effect_tables_read": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
