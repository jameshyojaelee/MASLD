#!/usr/bin/env bash
#SBATCH --job-name=model-check-329
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/model-check-328-21100337"
SOURCE_ARTIFACTS_SHA256=68d531468646c3f0a59a9193a1c82860de29b49a83aecbc7efbaa74b07d1ac1a
TARGET="${ROOT}/executions/model-check-329-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
RUNTIME="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738"
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON="${RUNTIME}/env/bin/python"

preserve_failure() {
  status=$?
  if [[ -d "${STAGE}" ]]; then
    if [[ "${status}" -ne 0 && ! -e "${TARGET}.failed" ]]; then
      mv "${STAGE}" "${TARGET}.failed"
    else
      rm -rf -- "${STAGE}"
    fi
  fi
  return "${status}"
}
trap preserve_failure EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Corgi+ readiness validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/runtime"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
export PYTHONPATH="${ROOT}/src"
cd "${ROOT}"
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_corgi_plus_current_terminal_readiness \
  > "${STAGE}/validation/unit-test.stdout.txt" 2>&1

"${PYTHON}" - "${SOURCE}" "${STAGE}/validation/validation.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

source = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(source)
audit = json.loads((source / "audit/audit.json").read_text())
assert audit["audit_status"] == "passed_terminal_blocked"
assert audit["current_primary_sources"]["corgi"]["high_level_corgi_plus_methods_are_placeholders"]
assert audit["current_primary_sources"]["corgi"]["tags"] == []
assert audit["current_primary_sources"]["corgi"]["current_head"] == "d84c00082da3d04eb05169e290e53f1c7da2b9c5"
assert audit["current_primary_sources"]["reproduction"]["current_head"] == "2367c42cf36630b7cf561eccc2e781abc22b36be"
assert audit["checkpoint"]["sha256"] == "01855614bfeffa72d5bcaaa28729c022ccfbc8f5fdcc0add8e04ba4533e94ccd"
assert audit["checkpoint"]["auxiliary_width"] == 26
assert audit["licenses"]["code_conflict"] == "UNRESOLVED"
assert not audit["licenses"]["open_champion_eligible"]
assert not audit["native_input_contract"]["native_input_bundle_released"]
assert not audit["native_input_contract"]["gene_count_matrix_is_valid_local_RNA"]
assert audit["reference_contract"]["biological_execution_gate"] == "blocked"
assert not any(row["corgi_plus_family_native_input_complete"] for row in audit["project_fixture_fit"].values())
assert audit["nearest_executable_lane"]["model_id"] == "corgi_regular"
assert audit["nearest_executable_lane"]["outcome_free_GSE296875_prediction_smoke_complete"]
assert not audit["nearest_executable_lane"]["relabel_as_Corgi_plus"]
assert not audit["context_borzoi_boundary"]["nearest_executable_now"]
assert not audit["terminal_disposition"]["production_or_GPU_item_authorized"]
assert not any(
    path.suffix.lower() in {".pt", ".pth", ".ckpt", ".safetensors", ".h5", ".hdf5"}
    for path in source.rglob("*")
    if path.is_file()
)
for field in (
    "checkpoint_deserialized_by_this_audit",
    "model_forward_executed_by_this_audit",
    "biological_data_read",
    "development_outcomes_read",
    "sealed_features_read",
    "sealed_outcomes_read",
    "gpu_job_submitted",
):
    assert audit[field] is False
result = {
    "schema_version": "masld-bench-corgi-plus-current-terminal-readiness-validation-v1",
    "source_artifact": source.name,
    "source_artifacts_sha256": "68d531468646c3f0a59a9193a1c82860de29b49a83aecbc7efbaa74b07d1ac1a",
    "unit_tests": 5,
    "frozen_tree_verified": True,
    "checkpoint_payload_absent": True,
    "biological_data_read": False,
    "development_outcomes_read": False,
    "sealed_features_read": False,
    "sealed_outcomes_read": False,
    "gpu_job_submitted": False,
    "terminal_disposition_reproduced": audit["terminal_disposition"],
    "nearest_executable_lane": audit["nearest_executable_lane"],
    "status": "passed",
}
Path(sys.argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
PY

"${PYTHON}" --version > "${STAGE}/runtime/python-version.txt" 2>&1
sha256sum \
  "${ROOT}/scripts/audit_corgi_plus_current_terminal_readiness.py" \
  "${ROOT}/tests/unit/test_audit_corgi_plus_current_terminal_readiness.py" \
  "${ROOT}/slurm/validate_corgi_plus_current_terminal_readiness_cpu.sbatch" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
validation = json.loads((root / "validation/validation.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "corgi_plus_current_terminal_readiness_validation",
        "model": "corgi_plus",
        "source_artifact": validation["source_artifact"],
        "source_artifacts_sha256": validation["source_artifacts_sha256"],
        "checkpoint_payload_absent": True,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_outcomes_read": False,
        "gpu_job_submitted": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Current Corgi+ terminal-readiness validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
