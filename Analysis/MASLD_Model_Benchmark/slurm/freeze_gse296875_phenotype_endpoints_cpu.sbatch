#!/usr/bin/env bash
#SBATCH --job-name=model-check-404
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
SUPPLEMENT=${ROOT}/executions/gse296875-donor-supplement-20260825-v2/supplement_lock
SUPPLEMENT_SHA256=cbc735fd1b0576c64ce0450d2b707b9dcc558166076f756df7505778136da8ee
CELLS=${SUPPLEMENT}/inventory/s1_donorinfo_cells.json
CELLS_SHA256=6e431a7b058376cf42adddcf0f1ff67c82545e5a666a15bd32f88007fb1f2cc8
JOIN_LOCK=${ROOT}/executions/gse296875-barcode-donor-join-20260825/join_lock
JOIN_LOCK_SHA256=f38e613f9263d24e51a01c678ddc7f8f3ca8781bb829ba08775f2c44dff9098f
TASKSPEC=${ROOT}/config/evaluation/gse296875_histopathology.toml
TASKSPEC_SHA256=1de972f98ce0ec67ad9d79612b735bee7d901a04066f7f5720b0cb26c24ad86d
DATASET=${ROOT}/config/datasets/gse296875.toml
DATASET_SHA256=9465b012bbf46611641555a5acae4d1f033d053c900bf944dcec02979c735b6f
TARGET=${ROOT}/executions/gse296875-phenotype-endpoints-20260825
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE296875 phenotype-endpoint artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${SUPPLEMENT_SHA256}" "${SUPPLEMENT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${CELLS_SHA256}" "${CELLS}" | sha256sum --check --strict
printf '%s  %s\n' "${JOIN_LOCK_SHA256}" "${JOIN_LOCK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${TASKSPEC_SHA256}" "${TASKSPEC}" | sha256sum --check --strict
printf '%s  %s\n' "${DATASET_SHA256}" "${DATASET}" | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse296875_phenotype_endpoints \
  tests.contract.test_gse296875_phenotype \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" scripts/freeze_gse296875_phenotype_endpoints.py \
  --cells "${CELLS}" \
  --join-lock "${JOIN_LOCK}" \
  --taskspec "${TASKSPEC}" \
  --output "${STAGE}/endpoint_lock" \
  > "${STAGE}/validation/endpoint_stdout.json" \
  2> "${STAGE}/validation/endpoint_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  config/evaluation/gse296875_histopathology.toml \
  config/datasets/gse296875.toml \
  scripts/freeze_gse296875_phenotype_endpoints.py \
  tests/unit/test_gse296875_phenotype_endpoints.py \
  slurm/freeze_gse296875_phenotype_endpoints_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
import tomllib

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root / "endpoint_lock")
receipt = json.loads((root / "endpoint_lock/endpoint_receipt.json").read_text())
rows = json.loads((root / "endpoint_lock/donor_endpoints.json").read_text())
with Path("config/evaluation/gse296875_histopathology.toml").open("rb") as handle:
    spec = tomllib.load(handle)
for row in rows:
    if not row["steatosis_observed"] and row["steatosis_numeric"] is not None:
        raise SystemExit("missing steatosis carries a value")
    if not row["fibrosis_observed"] and row["fibrosis_any"] is not None:
        raise SystemExit("missing fibrosis carries a value")
if (
    len(rows) != spec["analyzed_donors"]
    or receipt["steatosis_observed_donors"] != spec["steatosis"]["observed_donors"]
    or receipt["fibrosis_observed_donors"] != spec["fibrosis"]["observed_donors"]
    or receipt["missingness"] != "explicit_endpoint_mask_never_absence"
    or receipt["recoding"] != "none"
    or receipt["ordinal_stage_inference"] is not False
    or receipt["donors_under_18"] != spec["donors_under_18"]
    or manifest["metadata"]["champion_eligible"] is not False
):
    raise SystemExit("frozen endpoints differ from the TaskSpec")
freeze_tree(
    root,
    {
        "artifact_class": "gse296875_phenotype_endpoint_validation",
        "dataset_id": "gse296875",
        "analyzed_donors": len(rows),
        "steatosis_observed_donors": receipt["steatosis_observed_donors"],
        "fibrosis_observed_donors": receipt["fibrosis_observed_donors"],
        "champion_eligible": False,
        "status": "endpoints_and_masks_frozen",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/endpoint_lock/ARTIFACTS.json"
