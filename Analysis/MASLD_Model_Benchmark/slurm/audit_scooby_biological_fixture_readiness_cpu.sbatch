#!/usr/bin/env bash
#SBATCH --job-name=model-check-326
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/model-check-326-${SLURM_JOB_ID}"
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
  printf 'Refusing to overwrite Scooby biological-readiness audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/source-cache" "${STAGE}/validation" "${STAGE}/runtime"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
curl_flags=(--fail --location --silent --show-error --retry 6 --retry-all-errors --connect-timeout 30)

curl "${curl_flags[@]}" --output "${STAGE}/source-cache/scooby_onek1k.pinned.json" \
  https://huggingface.co/api/models/lauradmartens/onek1k-scooby/revision/dd7a2b7be54c906b43de4bc0159b169caaf057c5
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/scooby_onek1k.current.json" \
  https://huggingface.co/api/models/lauradmartens/onek1k-scooby
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/scooby_epicardioids.pinned.json" \
  https://huggingface.co/api/models/lauradmartens/epicardioids-scooby/revision/cefc3a626f5080baaf83878754ea924df6bd7e8d
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/scooby_epicardioids.current.json" \
  https://huggingface.co/api/models/lauradmartens/epicardioids-scooby
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/scooby_neurips.pinned.json" \
  https://huggingface.co/api/models/johahi/neurips-scooby/revision/1591d15696c65f14053e4495e52c74812a86fd8d
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/scooby_neurips.current.json" \
  https://huggingface.co/api/models/johahi/neurips-scooby
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/borzoi_converter.current.json" \
  https://huggingface.co/api/models/johahi/borzoi-replicate-0
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/scooby.LICENSE" \
  https://raw.githubusercontent.com/gagneurlab/scooby/dbfe6cd150d018d46a6d24ffa0406d7f2e89b471/LICENSE

export PYTHONPATH="${ROOT}/src"
cd "${ROOT}"
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_scooby_biological_fixture_readiness \
  > "${STAGE}/validation/unit-test.stdout.txt" 2>&1
"${PYTHON}" "${ROOT}/scripts/audit_scooby_biological_fixture_readiness.py" \
  --project-root "${ROOT}" \
  --source-cache "${STAGE}/source-cache" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json"

"${PYTHON}" - "${STAGE}/audit/audit.json" <<'PY'
from pathlib import Path
import json
import sys

audit = json.loads(Path(sys.argv[1]).read_text())
assert audit["audit_status"] == "passed_fail_closed"
assert audit["exact_code"]["license"] == "MIT"
assert audit["terms"]["canonical_calico_borzoi_weight_terms"] == "UNDECLARED"
assert audit["terms"]["merged_derivative_redistribution"] == "UNRESOLVED"
assert not audit["terms"]["open_champion_eligible"]
assert audit["dataset_fit"]["gse296875"]["supports_new_fold_fitted_context_model_in_principle"]
assert not audit["dataset_fit"]["gse296875"]["supports_frozen_released_scooby_context"]
assert audit["dataset_fit"]["gse281367"]["rna_state"] == "structurally_missing"
assert audit["central_dispatch"]["existing_item"].endswith("062-model-probe-284.json")
assert not audit["central_dispatch"]["duplicate_item_created"]
assert not audit["central_dispatch"]["manual_gpu_submission"]
assert audit["central_dispatch"]["qos"] == "nslab"
assert audit["central_dispatch"]["job_name"] == "model-probe-284"
assert not audit["checkpoint_deserialized"]
assert not audit["model_forward_executed"]
assert not audit["biological_matrices_read"]
assert not audit["project_outcomes_read"]
assert not audit["sealed_features_read"]
assert not audit["sealed_outcomes_read"]
assert not audit["gpu_job_submitted"]
PY

"${PYTHON}" --version > "${STAGE}/runtime/python-version.txt" 2>&1
"${PYTHON}" -m pip freeze --all > "${STAGE}/runtime/pip-freeze.txt"
uname -a > "${STAGE}/runtime/uname.txt"
env | LC_ALL=C sort | grep '^SLURM_' > "${STAGE}/runtime/slurm-environment.txt"
sha256sum \
  "${ROOT}/scripts/audit_scooby_biological_fixture_readiness.py" \
  "${ROOT}/tests/unit/test_audit_scooby_biological_fixture_readiness.py" \
  "${ROOT}/slurm/audit_scooby_biological_fixture_readiness_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "audit/audit.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "scooby_biological_fixture_readiness",
        "models": ["scooby_onek1k", "scooby_epicardioids", "scooby_neurips"],
        "audit_status": audit["audit_status"],
        "released_biological_smoke_authorized": False,
        "existing_synthetic_queue_item_verified": True,
        "duplicate_gpu_item_created": False,
        "gpu_job_submitted": False,
        "biological_matrices_read": False,
        "project_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_outcomes_read": False,
        "terminal_disposition": audit["terminal_disposition"],
        "status": "passed_fail_closed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Scooby biological-fixture readiness: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
