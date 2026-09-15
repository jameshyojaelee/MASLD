#!/usr/bin/env bash
#SBATCH --job-name=data-stage-403
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TASKSPEC=${ROOT}/config/evaluation/gse296875_histopathology.toml
TASKSPEC_SHA256=1de972f98ce0ec67ad9d79612b735bee7d901a04066f7f5720b0cb26c24ad86d
DATASET=${ROOT}/config/datasets/gse296875.toml
DATASET_SHA256=9465b012bbf46611641555a5acae4d1f033d053c900bf944dcec02979c735b6f
URL=https://ars.els-cdn.com/content/image/1-s2.0-S0002929725004343-mmc2.xlsx
WORKBOOK_NAME=1-s2.0-S0002929725004343-mmc2.xlsx
EXPECTED_SHA256=8356c2002059073faca89831181c8dc60d87857ebd408fe1021a7fabc15f73eb
EXPECTED_BYTES=5884547
TARGET=${ROOT}/executions/gse296875-donor-supplement-20260825-v2
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE296875 donor-supplement artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/download"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${TASKSPEC_SHA256}" "${TASKSPEC}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${DATASET_SHA256}" "${DATASET}" \
  | sha256sum --check --strict

# The publisher CDN is fetched with wget; the cluster has no curl egress.
# Some CDN edges refuse the default wget agent, so a descriptive agent is the
# second attempt.  Neither attempt can relax the frozen identity below.
DOWNLOAD=${STAGE}/download/${WORKBOOK_NAME}
AGENTS=("" "Mozilla/5.0 (X11; Linux x86_64) masld-bench-supplement-acquisition")
download_size_matched=false
for attempt in 0 1; do
  agent=${AGENTS[${attempt}]}
  wget_arguments=(
    --no-verbose
    --tries=5
    --waitretry=10
    --timeout=120
    --output-document="${DOWNLOAD}"
  )
  if [[ -n "${agent}" ]]; then
    wget_arguments+=(--user-agent="${agent}")
  fi
  if wget "${wget_arguments[@]}" "${URL}" \
    >> "${STAGE}/validation/wget_stdout.txt" \
    2>> "${STAGE}/validation/wget_stderr.txt"; then
    if [[ "$(stat -c '%s' "${DOWNLOAD}")" -eq "${EXPECTED_BYTES}" ]]; then
      download_size_matched=true
      break
    fi
  fi
  printf 'Attempt %s did not reach the frozen byte count.\n' "${attempt}" \
    >> "${STAGE}/validation/wget_stderr.txt"
done
if [[ "${download_size_matched}" != true ]]; then
  printf 'Supplement byte count differs: observed %s, frozen %s\n' \
    "$(stat -c '%s' "${DOWNLOAD}" 2>/dev/null || printf 0)" "${EXPECTED_BYTES}" >&2
  exit 1
fi
printf '%s  %s\n' "${EXPECTED_SHA256}" "${DOWNLOAD}" \
  | sha256sum --check --strict
sha256sum "${DOWNLOAD}" > "${STAGE}/validation/download.sha256"

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse296875_donor_supplement \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/acquire_gse296875_donor_supplement.py \
  --workbook "${DOWNLOAD}" \
  --taskspec "${TASKSPEC}" \
  --output "${STAGE}/supplement_lock" \
  > "${STAGE}/validation/acquisition_stdout.json" \
  2> "${STAGE}/validation/acquisition_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
wget --version | head -1 > "${STAGE}/validation/wget.version.txt"
sha256sum \
  config/evaluation/gse296875_histopathology.toml \
  config/datasets/gse296875.toml \
  scripts/acquire_gse296875_donor_supplement.py \
  tests/unit/test_gse296875_donor_supplement.py \
  slurm/acquire_gse296875_donor_supplement_io.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${EXPECTED_SHA256}" "${EXPECTED_BYTES}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
expected_sha256, expected_bytes = sys.argv[2], int(sys.argv[3])
manifest = verify_frozen_tree(root / "supplement_lock")
receipt = json.loads(
    (root / "supplement_lock/inventory/acquisition_receipt.json").read_text()
)
if (
    receipt["supplement_sha256"] != expected_sha256
    or receipt["supplement_size_bytes"] != expected_bytes
    or receipt["source_sheet"] != "S1 donorInfo"
    or receipt["recoding"] != "none"
    or manifest["metadata"]["champion_eligible"] is not False
):
    raise SystemExit("frozen supplement identity differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse296875_donor_supplement_acquisition_validation",
        "dataset_id": "gse296875",
        "supplement_sha256": expected_sha256,
        "supplement_size_bytes": expected_bytes,
        "champion_eligible": False,
        "external_or_sealed": False,
        "status": "acquired_hash_and_size_verified",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/supplement_lock/ARTIFACTS.json"
