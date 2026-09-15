#!/usr/bin/env bash
#SBATCH --job-name=transport-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/gse281367-atac-transport-production-readiness-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite transport audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse281367_atac_transport_production_readiness
PYTHONPATH="src:." python \
  scripts/audit_gse281367_atac_transport_production_readiness.py \
  --contract config/evaluation/gse281367_atac_transport_production_readiness.toml \
  --output "${STAGE}"
sha256sum \
  config/evaluation/gse281367_atac_transport_production_readiness.toml \
  scripts/audit_gse281367_atac_transport_production_readiness.py \
  tests/unit/test_audit_gse281367_atac_transport_production_readiness.py \
  config/evaluation/gse281367_activation_readiness.toml \
  config/evaluation/observed_multiome_task.toml \
  config/evaluation/observed_multiome_development_gate.toml \
  OVERALL_PLAN.md \
  slurm/audit_gse281367_atac_transport_production_readiness_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src:." python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse281367_atac_transport_production_readiness_audit",
            "outcome_visibility": "metadata_and_window_coordinates_only",
            "sealed_outcomes_read": False,
            "biological_outcome_arrays_read": False,
            "condition_values_read": False,
            "frozen_outcome_authority_valid": True,
            "label_agnostic_production_ready": False,
            "observed_atac_production_ready": False,
            "role": "project_exposed_atac_only_development",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE281367 ATAC transport production-readiness audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
