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
TARGET="${ROOT}/executions/gse281367-activation-readiness-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite cohort audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse281367_activation_readiness
PYTHONPATH="src:." python scripts/audit_gse281367_activation_readiness.py \
  --contract config/evaluation/gse281367_activation_readiness.toml \
  --output "${STAGE}"
sha256sum \
  config/evaluation/gse281367_activation_readiness.toml \
  config/datasets/gse281367.toml \
  scripts/audit_gse281367_activation_readiness.py \
  tests/unit/test_audit_gse281367_activation_readiness.py \
  config/evaluation/cross_cohort_expansion.toml \
  config/evaluation/observed_multiome_task.toml \
  OVERALL_PLAN.md \
  slurm/audit_gse281367_activation_readiness_cpu.sbatch \
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
            "artifact_class": "gse281367_activation_readiness_audit",
            "outcome_visibility": "public_metadata_and_frozen_audit_only",
            "sealed_outcomes_read": False,
            "biological_matrices_downloaded": False,
            "controlled_data_accessed": False,
            "activates_source_native_model_input": False,
            "fixed_window_transport_active": True,
            "role": "project_exposed_atac_only_development",
            "source_label_semantics": "conflict_geo_mash_article_masld",
            "rna_state": "structurally_missing",
            "external_or_sealed_champion_eligible": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE281367 activation-readiness audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
