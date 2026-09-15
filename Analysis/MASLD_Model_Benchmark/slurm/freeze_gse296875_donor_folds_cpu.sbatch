#!/usr/bin/env bash
#SBATCH --job-name=model-check-405
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
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
ENDPOINT_LOCK=${ROOT}/executions/gse296875-phenotype-endpoints-20260825/endpoint_lock
ENDPOINT_LOCK_SHA256=6e8c1bde8e85d2ff5369aee43141d0f524897867df4440a9330ebe89dfb6b477
JOIN_LOCK=${ROOT}/executions/gse296875-barcode-donor-join-20260825/join_lock
JOIN_LOCK_SHA256=f38e613f9263d24e51a01c678ddc7f8f3ca8781bb829ba08775f2c44dff9098f
MEMBERSHIP=${ROOT}/executions/gse296875-fragment-membership-21063829
MEMBERSHIP_SHA256=582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64
MEMBERSHIP_FOLDS_SHA256=90b21fa6418baa344722526d580a4a4028775b50ff1374ffc6829c99f9ed36d4
TASKSPEC=${ROOT}/config/evaluation/gse296875_histopathology.toml
TASKSPEC_SHA256=1de972f98ce0ec67ad9d79612b735bee7d901a04066f7f5720b0cb26c24ad86d
DATASET=${ROOT}/config/datasets/gse296875.toml
DATASET_SHA256=9465b012bbf46611641555a5acae4d1f033d053c900bf944dcec02979c735b6f
TARGET=${ROOT}/executions/gse296875-donor-folds-20260825
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE296875 donor-fold artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${ENDPOINT_LOCK_SHA256}" "${ENDPOINT_LOCK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${JOIN_LOCK_SHA256}" "${JOIN_LOCK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${MEMBERSHIP_SHA256}" "${MEMBERSHIP}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${MEMBERSHIP_FOLDS_SHA256}" "${MEMBERSHIP}/donor_folds.tsv" | sha256sum --check --strict
printf '%s  %s\n' "${TASKSPEC_SHA256}" "${TASKSPEC}" | sha256sum --check --strict
printf '%s  %s\n' "${DATASET_SHA256}" "${DATASET}" | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse296875_donor_folds \
  tests.contract.test_gse296875_phenotype \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/freeze_gse296875_donor_folds.py \
  --endpoint-lock "${ENDPOINT_LOCK}" \
  --join-lock "${JOIN_LOCK}" \
  --membership "${MEMBERSHIP}" \
  --taskspec "${TASKSPEC}" \
  --outer-folds 5 \
  --output "${STAGE}/split_lock" \
  > "${STAGE}/validation/split_stdout.json" \
  2> "${STAGE}/validation/split_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  config/evaluation/gse296875_histopathology.toml \
  config/datasets/gse296875.toml \
  scripts/freeze_gse296875_donor_folds.py \
  tests/unit/test_gse296875_donor_folds.py \
  slurm/freeze_gse296875_donor_folds_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys
import tomllib

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root / "split_lock")
receipt = json.loads((root / "split_lock/split_receipt.json").read_text())
with Path("config/evaluation/gse296875_histopathology.toml").open("rb") as handle:
    spec = tomllib.load(handle)
with (root / "split_lock/donor_folds.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
donors = [row["donor_id"] for row in rows]
if len(set(donors)) != len(donors) or len(donors) != spec["analyzed_donors"]:
    raise SystemExit("donor-grouped split does not cover each donor exactly once")
if sum(1 for row in rows if row["in_adult_only_sensitivity"] == "false") != spec["donors_under_18"]:
    raise SystemExit("the adult-only mask does not match the frozen under-18 count")
if (
    receipt["wells"] != spec["technical_batches"]
    or receipt["frozen_before_any_fitting"] is not True
    or receipt["under_18_donors_are_masked_not_dropped"] is not True
    or receipt["donor_is_nested_within_well"] is not True
    or manifest["metadata"]["champion_eligible"] is not False
):
    raise SystemExit("frozen split differs from the TaskSpec")
freeze_tree(
    root,
    {
        "artifact_class": "gse296875_donor_grouped_split_validation",
        "dataset_id": "gse296875",
        "donors": receipt["donors"],
        "wells": receipt["wells"],
        "outer_folds": receipt["outer_folds"],
        "pooled_out_of_fold_metric_required": receipt["pooled_out_of_fold_metric_required"],
        "champion_eligible": False,
        "status": "splits_frozen_before_fitting",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/split_lock/ARTIFACTS.json"
