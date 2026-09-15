#!/usr/bin/env bash
#SBATCH --job-name=model-data-402
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
CONTINUAL=${PROJECT_ROOT}/Analysis/SingleCell/continual_integration
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
JOIN=${CONTINUAL}/candidates/external-gse296875-metadata-v45-2026-08-17T230000Z/mapping_metadata.tsv
JOIN_SHA256=5d3d946607f55b140ffa4a44369cf35f6e9e91f0f17276afc0e3e9f008554781
SOURCE_LOCK=${CONTINUAL}/reference/gse296875_source_lock_v44.json
SOURCE_LOCK_SHA256=6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620
RAW_RNA=${CONTINUAL}/reference/gse296875_raw_rna_sources_v44.tsv
RAW_RNA_SHA256=610f4ec7f71abf5d9e049b1093f69bfad85c02ebf8570e44882890fac15941ae
MEMBERSHIP=${ROOT}/executions/gse296875-fragment-membership-21063829
MEMBERSHIP_SHA256=582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64
DATASET=${ROOT}/config/datasets/gse296875.toml
DATASET_SHA256=9465b012bbf46611641555a5acae4d1f033d053c900bf944dcec02979c735b6f
TASKSPEC=${ROOT}/config/evaluation/gse296875_histopathology.toml
TASKSPEC_SHA256=1de972f98ce0ec67ad9d79612b735bee7d901a04066f7f5720b0cb26c24ad86d
TARGET=${ROOT}/executions/gse296875-barcode-donor-join-20260825
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE296875 barcode-join artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${JOIN_SHA256}" "${JOIN}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_LOCK_SHA256}" "${SOURCE_LOCK}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RAW_RNA_SHA256}" "${RAW_RNA}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${MEMBERSHIP_SHA256}" "${MEMBERSHIP}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${DATASET_SHA256}" "${DATASET}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${TASKSPEC_SHA256}" "${TASKSPEC}" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse296875_barcode_donor_join \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/freeze_gse296875_barcode_donor_join.py \
  --join "${JOIN}" \
  --source-lock "${SOURCE_LOCK}" \
  --raw-rna-sources "${RAW_RNA}" \
  --membership "${MEMBERSHIP}" \
  --expected-donors 39 \
  --expected-wells 8 \
  --output "${STAGE}/join_lock" \
  > "${STAGE}/validation/join_stdout.json" \
  2> "${STAGE}/validation/join_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  config/datasets/gse296875.toml \
  config/evaluation/gse296875_histopathology.toml \
  scripts/freeze_gse296875_barcode_donor_join.py \
  tests/unit/test_gse296875_barcode_donor_join.py \
  slurm/freeze_gse296875_barcode_donor_join_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${JOIN_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
import tomllib

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
join_sha256 = sys.argv[2]
manifest = verify_frozen_tree(root / "join_lock")
receipt = json.loads((root / "join_lock/join_receipt.json").read_text())
with Path("config/datasets/gse296875.toml").open("rb") as handle:
    dataset = tomllib.load(handle)
if (
    receipt["join_sha256"] != join_sha256
    or receipt["donors"] != dataset["expected_biological_units"]
    or receipt["wells"] != 8
    or receipt["nuclei"] != 68398
    or receipt["unique_cell_id"] != receipt["nuclei"]
    or receipt["unique_well_namespaced_barcode"] != receipt["nuclei"]
    or receipt["biological_unit"] != "donor"
    or receipt["technical_batch"] != "well"
    or manifest["metadata"]["champion_eligible"] is not False
):
    raise SystemExit("frozen barcode join differs from the dataset contract")
freeze_tree(
    root,
    {
        "artifact_class": "gse296875_barcode_donor_join_validation",
        "dataset_id": "gse296875",
        "join_sha256": join_sha256,
        "donors": receipt["donors"],
        "wells": receipt["wells"],
        "nuclei": receipt["nuclei"],
        "donor_is_nested_within_well": receipt["donor_is_nested_within_well"],
        "well_namespacing_is_load_bearing": receipt["well_namespacing_is_load_bearing"],
        "champion_eligible": False,
        "status": "join_frozen_and_hashed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/join_lock/ARTIFACTS.json"
