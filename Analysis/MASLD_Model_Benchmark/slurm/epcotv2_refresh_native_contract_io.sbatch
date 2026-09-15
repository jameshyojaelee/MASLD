#!/usr/bin/env bash
#SBATCH --job-name=model-audit-112
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/model-audit-112-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
ADMISSION="${ROOT}/executions/epibert-epcotv2-admission-sources-21064417"
LEGACY_FIXTURE="${ROOT}/executions/epibert-epcotv2-no-outcome-fixture-21064487"
SPACE_REVISION=f93a2a944a49797a111f0294709416d38bedbcbe
REFERENCE_REVISION=f445058917800811cb82c38b2c806674cc02b47f
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
  printf 'Refusing to overwrite EPCOTv2 native-contract artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/terms" "${STAGE}/validation"

curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30)
curl "${curl_flags[@]}" --output "${STAGE}/terms/space_current.json" \
  https://huggingface.co/api/spaces/luosanj/EPCOTv2
curl "${curl_flags[@]}" --output "${STAGE}/terms/space_pinned.json" \
  "https://huggingface.co/api/spaces/luosanj/EPCOTv2/revision/${SPACE_REVISION}"
curl "${curl_flags[@]}" --output "${STAGE}/terms/reference_current.json" \
  https://huggingface.co/api/datasets/luosanj/epcotv2_data
curl "${curl_flags[@]}" --output "${STAGE}/terms/reference_pinned.json" \
  "https://huggingface.co/api/datasets/luosanj/epcotv2_data/revision/${REFERENCE_REVISION}"

space_pinned_status=$(curl --location --silent --show-error --output /dev/null --write-out '%{http_code}' \
  "https://huggingface.co/spaces/luosanj/EPCOTv2/resolve/${SPACE_REVISION}/LICENSE")
space_current_status=$(curl --location --silent --show-error --output /dev/null --write-out '%{http_code}' \
  https://huggingface.co/spaces/luosanj/EPCOTv2/resolve/main/LICENSE)
reference_pinned_status=$(curl --location --silent --show-error --output /dev/null --write-out '%{http_code}' \
  "https://huggingface.co/datasets/luosanj/epcotv2_data/resolve/${REFERENCE_REVISION}/LICENSE")
reference_current_status=$(curl --location --silent --show-error --output /dev/null --write-out '%{http_code}' \
  https://huggingface.co/datasets/luosanj/epcotv2_data/resolve/main/LICENSE)
printf '{"reference_current":%s,"reference_pinned":%s,"space_current":%s,"space_pinned":%s}\n' \
  "${reference_current_status}" "${reference_pinned_status}" "${space_current_status}" "${space_pinned_status}" \
  > "${STAGE}/terms/license_endpoint_status.json"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
export PYTHONPATH="${ROOT}/src"
cd "${ROOT}"
"${PYTHON}" -m unittest -v tests.unit.test_audit_epcotv2_native_contract \
  > "${STAGE}/validation/unit-test.stdout.txt" 2>&1
"${PYTHON}" "${ROOT}/scripts/audit_epcotv2_native_contract.py" \
  --admission "${ADMISSION}" \
  --legacy-fixture "${LEGACY_FIXTURE}" \
  --terms-snapshot "${STAGE}/terms" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json"

"${PYTHON}" - "${STAGE}/audit/audit.json" <<'PY'
from pathlib import Path
import json
import sys

audit = json.loads(Path(sys.argv[1]).read_text())
assert audit["status"] == "pass"
assert audit["terms"]["terminal_disposition"] == "blocked_terms_no_checkpoint_download_or_execution"
assert audit["fixture_scope"]["source_defined_preprocessing_fixture_without_biological_channel_names"] == "admissible"
assert audit["fixture_scope"]["checkpoint_forward_or_numeric_fixture"].startswith("blocked_")
assert not audit["project_data_read"] and not audit["outcomes_read"] and not audit["sealed_labels_read"]
PY

sha256sum \
  "${ROOT}/scripts/audit_epcotv2_native_contract.py" \
  "${ROOT}/tests/unit/test_audit_epcotv2_native_contract.py" \
  "${ROOT}/slurm/epcotv2_refresh_native_contract_io.sbatch" \
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
        "artifact_class": "epcotv2_native_contract_audit",
        "model": "epcotv2",
        "checkpoint_downloaded": False,
        "model_forward_executed": False,
        "preprocessing_fixture_status": "passed_biological_free",
        "terminal_disposition": audit["terms"]["terminal_disposition"],
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_labels_read": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EPCOTv2 native-contract audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
