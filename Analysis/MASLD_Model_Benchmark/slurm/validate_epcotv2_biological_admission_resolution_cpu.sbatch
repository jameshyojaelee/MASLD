#!/usr/bin/env bash
#SBATCH --job-name=model-check-320
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
SOURCE="${ROOT}/executions/model-check-319-21100063"
SOURCE_ARTIFACTS_SHA256=4fe7ac3dedc49c807cc8d5e2f8a12444dbe88e29cd6cf83c3335acf4142e8c61
TARGET="${ROOT}/executions/model-check-320-${SLURM_JOB_ID}"
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
  printf 'Refusing to overwrite EPCOTv2 validation artifact: %s\n' "${TARGET}" >&2
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
  tests.unit.test_audit_epcotv2_biological_admission_resolution \
  > "${STAGE}/validation/unit-test.stdout.txt" 2>&1
"${PYTHON}" - "${SOURCE}" "${STAGE}/validation/validation.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

source = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(source)
audit = json.loads((source / "audit/audit.json").read_text())
assert audit["audit_status"] == "passed_fail_closed"
assert audit["publication_tutorial_sequence_contract"]["base_to_channel"] == {
    "A": 0, "C": 1, "G": 2, "T": 3
}
assert audit["publication_tutorial_sequence_contract"]["neutral_base"] == {
    "neutral_alphabet": ["N"], "neutral_value_per_channel": 0.25
}
assert audit["space_reference_array_contract"]["base_to_row_mapping"] == "UNRESOLVED"
assert audit["space_reference_array_contract"]["checkpoint_binding"] == "UNRESOLVED"
assert audit["reference_identity"]["exact_checkpoint_reference_identity"] == "UNRESOLVED"
assert audit["terms"]["checkpoint_weights"] == "UNDECLARED"
assert audit["terms"]["reference_arrays"] == "UNDECLARED"
assert audit["released_execution_path"]["safe_local_execution_path"] == "UNAVAILABLE"
assert audit["terminal_disposition"]["biological_execution"] == "blocked"
assert audit["terminal_disposition"]["champion_eligibility"] == "ineligible"
assert audit["terminal_disposition"]["gpu_admission"] is False
assert not audit["checkpoint_object_present_in_audit"]
assert not any(path.suffix in {".pt", ".pth", ".ckpt", ".safetensors"} for path in source.rglob("*"))
assert not audit["project_data_read"] and not audit["outcomes_read"]
assert not audit["sealed_labels_read"]
assert not audit["upstream_code_imported_or_executed"]
assert not audit["pickle_or_torch_deserialization_performed"]
result = {
    "schema_version": "masld-bench-epcotv2-biological-admission-validation-v1",
    "source_artifact": source.name,
    "source_artifacts_sha256": "4fe7ac3dedc49c807cc8d5e2f8a12444dbe88e29cd6cf83c3335acf4142e8c61",
    "unit_tests": 5,
    "frozen_tree_verified": True,
    "checkpoint_payload_absent": True,
    "outcomes_read": False,
    "sealed_labels_read": False,
    "terminal_disposition_reproduced": "blocked_biological_execution",
    "status": "passed",
}
Path(sys.argv[2]).write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
PY

"${PYTHON}" --version > "${STAGE}/runtime/python-version.txt" 2>&1
sha256sum \
  "${ROOT}/scripts/audit_epcotv2_biological_admission_resolution.py" \
  "${ROOT}/tests/unit/test_audit_epcotv2_biological_admission_resolution.py" \
  "${ROOT}/slurm/validate_epcotv2_biological_admission_resolution_cpu.sbatch" \
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
        "artifact_class": "epcotv2_biological_admission_validation",
        "model": "epcotv2",
        "source_artifact": validation["source_artifact"],
        "source_artifacts_sha256": validation["source_artifacts_sha256"],
        "checkpoint_payload_absent": True,
        "outcomes_read": False,
        "sealed_labels_read": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EPCOTv2 biological-admission validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
