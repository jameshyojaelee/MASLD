#!/usr/bin/env bash
#SBATCH --job-name=cohort-audit
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
TARGET="${ROOT}/executions/external-observed-multiome-cohort-audit-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite cohort audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_external_observed_multiome_cohorts
PYTHONPATH="src:." python scripts/audit_external_observed_multiome_cohorts.py \
  --contract config/evaluation/external_observed_multiome_cohort_expansion.toml \
  --output "${STAGE}"
sha256sum \
  config/evaluation/external_observed_multiome_cohort_expansion.toml \
  scripts/audit_external_observed_multiome_cohorts.py \
  tests/unit/test_audit_external_observed_multiome_cohorts.py \
  config/datasets/fnih_86_liver.toml \
  config/datasets/gse244832.toml \
  config/datasets/gse296875.toml \
  config/evaluation/observed_multiome_task.toml \
  slurm/audit_external_observed_multiome_cohorts_cpu.sbatch \
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
            "artifact_class": "external_observed_multiome_cohort_expansion_audit",
            "outcome_visibility": "public_metadata_only",
            "sealed_outcomes_read": False,
            "biological_matrices_downloaded": False,
            "activates_dataset": False,
            "anchor_excluded_from_ranking": "gse296875",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'External observed-multiome cohort audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
