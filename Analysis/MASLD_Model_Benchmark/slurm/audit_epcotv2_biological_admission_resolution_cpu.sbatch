#!/usr/bin/env bash
#SBATCH --job-name=model-check-319
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/model-check-319-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
ADMISSION="${ROOT}/executions/epibert-epcotv2-admission-sources-21064417"
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
  printf 'Refusing to overwrite EPCOTv2 resolution artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/source-cache/downloads" \
  "${STAGE}/source-cache/metadata" \
  "${STAGE}/source-cache/sources/space/curriculum/loralib" \
  "${STAGE}/source-cache/sources/space/erna" \
  "${STAGE}/source-cache/sources/reference" \
  "${STAGE}/validation" \
  "${STAGE}/runtime"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
curl_flags=(--fail --location --silent --show-error --retry 6 --retry-all-errors --connect-timeout 30)

curl "${curl_flags[@]}" --output "${STAGE}/source-cache/metadata/space-current.json" \
  https://huggingface.co/api/spaces/luosanj/EPCOTv2
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/metadata/space-pinned.json" \
  "https://huggingface.co/api/spaces/luosanj/EPCOTv2/revision/${SPACE_REVISION}"
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/metadata/reference-current.json" \
  https://huggingface.co/api/datasets/luosanj/epcotv2_data
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/metadata/reference-pinned.json" \
  "https://huggingface.co/api/datasets/luosanj/epcotv2_data/revision/${REFERENCE_REVISION}"
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/metadata/kipoiseq-0.5.2-pypi.json" \
  https://pypi.org/pypi/kipoiseq/0.5.2/json

"${PYTHON}" - "${STAGE}/source-cache/metadata/kipoiseq-0.5.2-pypi.json" \
  "${STAGE}/source-cache/metadata/kipoiseq-sdist-url.txt" <<'PY'
from pathlib import Path
import json
import sys

metadata = json.loads(Path(sys.argv[1]).read_text())
rows = [row for row in metadata["urls"] if row.get("packagetype") == "sdist"]
if len(rows) != 1:
    raise SystemExit("unexpected kipoiseq sdist cardinality")
Path(sys.argv[2]).write_text(rows[0]["url"] + "\n")
PY
read -r kipoiseq_sdist_url < "${STAGE}/source-cache/metadata/kipoiseq-sdist-url.txt"
curl "${curl_flags[@]}" --output "${STAGE}/source-cache/downloads/kipoiseq-0.5.2.tar.gz" \
  "${kipoiseq_sdist_url}"
printf '%s  %s\n' \
  b1d8f717326522e4fd88baca8bb288050a6da2ec75e23b82006c55c6bc68204d \
  "${STAGE}/source-cache/downloads/kipoiseq-0.5.2.tar.gz" \
  | sha256sum --check --strict

space_files=(
  README.md
  app.py
  curriculum/lora_prompt_model.py
  curriculum/loralib/__init__.py
  curriculum/loralib/layers.py
  curriculum/loralib/utils.py
  erna/util.py
  func_gradio.py
  requirements.txt
)
for relative in "${space_files[@]}"; do
  curl "${curl_flags[@]}" \
    --output "${STAGE}/source-cache/sources/space/${relative}" \
    "https://huggingface.co/spaces/luosanj/EPCOTv2/resolve/${SPACE_REVISION}/${relative}"
done
curl "${curl_flags[@]}" \
  --output "${STAGE}/source-cache/sources/reference/README.md" \
  "https://huggingface.co/datasets/luosanj/epcotv2_data/resolve/${REFERENCE_REVISION}/README.md"

space_pinned_status=$(curl --location --silent --show-error --output /dev/null --write-out '%{http_code}' \
  "https://huggingface.co/spaces/luosanj/EPCOTv2/resolve/${SPACE_REVISION}/LICENSE")
space_current_status=$(curl --location --silent --show-error --output /dev/null --write-out '%{http_code}' \
  https://huggingface.co/spaces/luosanj/EPCOTv2/resolve/main/LICENSE)
reference_pinned_status=$(curl --location --silent --show-error --output /dev/null --write-out '%{http_code}' \
  "https://huggingface.co/datasets/luosanj/epcotv2_data/resolve/${REFERENCE_REVISION}/LICENSE")
reference_current_status=$(curl --location --silent --show-error --output /dev/null --write-out '%{http_code}' \
  https://huggingface.co/datasets/luosanj/epcotv2_data/resolve/main/LICENSE)
printf '{"reference_current":%s,"reference_pinned":%s,"space_current":%s,"space_pinned":%s}\n' \
  "${reference_current_status}" "${reference_pinned_status}" \
  "${space_current_status}" "${space_pinned_status}" \
  > "${STAGE}/source-cache/metadata/license-endpoint-status.json"

export PYTHONPATH="${ROOT}/src"
cd "${ROOT}"
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_epcotv2_biological_admission_resolution \
  > "${STAGE}/validation/unit-test.stdout.txt" 2>&1
"${PYTHON}" "${ROOT}/scripts/audit_epcotv2_biological_admission_resolution.py" \
  --admission "${ADMISSION}" \
  --source-cache "${STAGE}/source-cache" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json"

"${PYTHON}" - "${STAGE}/audit/audit.json" <<'PY'
from pathlib import Path
import json
import sys

audit = json.loads(Path(sys.argv[1]).read_text())
assert audit["audit_status"] == "passed_fail_closed"
assert audit["publication_tutorial_sequence_contract"]["base_to_channel"] == {
    "A": 0, "C": 1, "G": 2, "T": 3
}
assert audit["space_reference_array_contract"]["base_to_row_mapping"] == "UNRESOLVED"
assert audit["reference_identity"]["current_external_checksum_used_for_admission"] is False
assert audit["terminal_disposition"]["biological_execution"] == "blocked"
assert audit["terminal_disposition"]["gpu_admission"] is False
assert not audit["checkpoint_object_present_in_audit"]
assert not audit["project_data_read"] and not audit["outcomes_read"]
assert not audit["sealed_labels_read"]
assert not audit["upstream_code_imported_or_executed"]
assert not audit["pickle_or_torch_deserialization_performed"]
PY

"${PYTHON}" --version > "${STAGE}/runtime/python-version.txt" 2>&1
"${PYTHON}" -m pip freeze --all > "${STAGE}/runtime/pip-freeze.txt"
uname -a > "${STAGE}/runtime/uname.txt"
env | LC_ALL=C sort | grep '^SLURM_' > "${STAGE}/runtime/slurm-environment.txt"
sha256sum \
  "${ROOT}/scripts/audit_epcotv2_biological_admission_resolution.py" \
  "${ROOT}/tests/unit/test_audit_epcotv2_biological_admission_resolution.py" \
  "${ROOT}/slurm/audit_epcotv2_biological_admission_resolution_cpu.sbatch" \
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
        "artifact_class": "epcotv2_biological_admission_resolution",
        "model": "epcotv2",
        "checkpoint_downloaded": False,
        "checkpoint_deserialized": False,
        "upstream_code_executed": False,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_labels_read": False,
        "terminal_disposition": audit["terminal_disposition"],
        "status": "passed_fail_closed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EPCOTv2 biological-admission resolution: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
