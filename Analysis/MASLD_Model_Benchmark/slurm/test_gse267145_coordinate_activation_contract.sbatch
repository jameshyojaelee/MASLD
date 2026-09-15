#!/usr/bin/env bash
#SBATCH --job-name=gse267145-contract
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
TARGET="${ROOT}/executions/gse267145-coordinate-contract-test-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 contract test: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_coordinate_activation \
  tests.contract.test_gse267145_coordinate_activation_contract \
  tests.contract.test_cross_cohort_plan \
  tests.contract.test_registry_files.RegistryFileContractTest.test_all_registry_files_parse_and_have_expected_schema
python -m pip freeze --all > "${STAGE}/pip-freeze.txt"
python - "${STAGE}/test_receipt.json" <<'PY'
from pathlib import Path
import json
import sys

Path(sys.argv[1]).write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-gse267145-coordinate-contract-test-v1",
            "status": "pass",
            "test_modules": [
                "tests.unit.test_audit_gse267145_coordinate_activation",
                "tests.contract.test_gse267145_coordinate_activation_contract",
                "tests.contract.test_cross_cohort_plan",
                "tests.contract.test_registry_files.RegistryFileContractTest.test_all_registry_files_parse_and_have_expected_schema",
            ],
            "activation_artifacts_sha256": "e5042eb43dddd84cda81c5aeed5d6dbeb0378ab1c942a0183e2ef1b41cd18b20",
            "coordinate_semantics_resolved": False,
            "sequence_extraction_allowed": False,
            "task_specific_coordinate_independent_training_allowed": True,
            "champion_eligible": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
PY
sha256sum \
  scripts/audit_gse267145_coordinate_activation.py \
  tests/unit/test_audit_gse267145_coordinate_activation.py \
  tests/contract/test_gse267145_coordinate_activation_contract.py \
  tests/contract/test_cross_cohort_plan.py \
  tests/contract/test_registry_files.py \
  slurm/test_gse267145_coordinate_activation_contract.sbatch \
  config/datasets/gse267145_znf469_human_liver.toml \
  config/evaluation/cross_cohort_expansion.toml \
  config/evaluation/family_native_tournament.toml \
  config/evaluation/paired_bulk_rna_h3k27ac_task.toml \
  config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json \
  OVERALL_PLAN.md \
  > "${STAGE}/source.sha256"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]).resolve(strict=True),
    {
        "artifact_class": "gse267145_coordinate_activation_contract_test",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "task_id": "paired_bulk_rna_h3k27ac",
        "activation_artifacts_sha256": "e5042eb43dddd84cda81c5aeed5d6dbeb0378ab1c942a0183e2ef1b41cd18b20",
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 coordinate activation contract tests: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
