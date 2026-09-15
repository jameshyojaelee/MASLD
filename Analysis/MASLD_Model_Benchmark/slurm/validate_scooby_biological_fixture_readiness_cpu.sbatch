#!/usr/bin/env bash
#SBATCH --job-name=model-check-327
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
SOURCE="${ROOT}/executions/model-check-326-21100204"
SOURCE_ARTIFACTS_SHA256=d9ffa5151a5cb8dc91d69fcf7580497b237cc3cdb5c9adb3a98a7f9f9aa51b68
TARGET="${ROOT}/executions/model-check-327-${SLURM_JOB_ID}"
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
  printf 'Refusing to overwrite Scooby readiness validation: %s\n' "${TARGET}" >&2
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
  tests.unit.test_audit_scooby_biological_fixture_readiness \
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
assert set(audit["released_checkpoints"]) == {
    "scooby_onek1k", "scooby_epicardioids", "scooby_neurips"
}
assert audit["released_checkpoints"]["scooby_epicardioids"]["context_dim"] == 50
assert audit["released_checkpoints"]["scooby_neurips"]["context_dim"] == 14
assert audit["released_checkpoints"]["scooby_onek1k"]["context_dim"] == 10
assert audit["released_checkpoints"]["scooby_neurips"]["pinned_terms"]["terms_declared"] is False
assert audit["terms"]["borzoi_converter_current_snapshot"]["declared_license"] == "cc-by-4.0"
assert audit["terms"]["canonical_calico_borzoi_weight_terms"] == "UNDECLARED"
assert audit["native_contract"]["exact_native_FASTA_identity"] == "UNRESOLVED"
assert audit["native_contract"]["reference_and_v32_to_v49_biological_execution_gate"] == "blocked"
assert audit["biological_fixture_decisions"]["released_checkpoint_on_gse296875"]["status"] == "blocked"
assert audit["biological_fixture_decisions"]["project_fitted_scooby_adaptation_on_gse296875"]["status"] == "conditional_not_authorized"
assert audit["biological_fixture_decisions"]["released_checkpoint_on_gse281367"]["status"] == "not_applicable"
assert audit["central_dispatch"]["enabled"]
assert not audit["central_dispatch"]["duplicate_item_created"]
assert not audit["central_dispatch"]["manual_gpu_submission"]
assert audit["central_dispatch"]["qos"] == "nslab"
assert not any(
    path.suffix.lower() in {".pt", ".pth", ".ckpt", ".safetensors", ".h5", ".hdf5"}
    for path in source.rglob("*")
    if path.is_file()
)
for field in (
    "checkpoint_deserialized",
    "model_forward_executed",
    "biological_matrices_read",
    "project_outcomes_read",
    "sealed_features_read",
    "sealed_outcomes_read",
    "gpu_job_submitted",
):
    assert audit[field] is False
result = {
    "schema_version": "masld-bench-scooby-biological-fixture-readiness-validation-v1",
    "source_artifact": source.name,
    "source_artifacts_sha256": "d9ffa5151a5cb8dc91d69fcf7580497b237cc3cdb5c9adb3a98a7f9f9aa51b68",
    "unit_tests": 5,
    "frozen_tree_verified": True,
    "checkpoint_payload_absent": True,
    "biological_matrices_read": False,
    "project_outcomes_read": False,
    "sealed_features_read": False,
    "sealed_outcomes_read": False,
    "gpu_job_submitted": False,
    "terminal_disposition_reproduced": audit["terminal_disposition"],
    "status": "passed",
}
Path(sys.argv[2]).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
PY

"${PYTHON}" --version > "${STAGE}/runtime/python-version.txt" 2>&1
sha256sum \
  "${ROOT}/scripts/audit_scooby_biological_fixture_readiness.py" \
  "${ROOT}/tests/unit/test_audit_scooby_biological_fixture_readiness.py" \
  "${ROOT}/slurm/validate_scooby_biological_fixture_readiness_cpu.sbatch" \
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
        "artifact_class": "scooby_biological_fixture_readiness_validation",
        "models": ["scooby_onek1k", "scooby_epicardioids", "scooby_neurips"],
        "source_artifact": validation["source_artifact"],
        "source_artifacts_sha256": validation["source_artifacts_sha256"],
        "checkpoint_payload_absent": True,
        "biological_matrices_read": False,
        "project_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_outcomes_read": False,
        "gpu_job_submitted": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Scooby biological-fixture readiness validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
