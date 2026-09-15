#!/usr/bin/env bash
#SBATCH --job-name=data-check-248
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-check-248-%j.out
#SBATCH --error=executions/logs/data-check-248-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
REFERENCE=${ROOT}/executions/gse105127-production-9675cd48993468d4/reference
REFERENCE_SHA256=9575aa24dc39f855b961eecf53e29243b4a45fca6770c7337e05901bc7d1ec89
REFERENCE_RECEIPT_SHA256=5a518bc6bf4ae84e2a320b810cf85ace8fa517cd91200194d728d2cadeadddb6
TARGET=${ROOT}/executions/data-check-248-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE105127 reference-status validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${REFERENCE_SHA256}" "${REFERENCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${REFERENCE_RECEIPT_SHA256}" "${REFERENCE}/receipt.json" \
  | sha256sum --check --strict

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse105127_reference_status_repair \
  > "${STAGE}/unittest.stdout.txt" \
  2> "${STAGE}/unittest.stderr.txt"
"${PYTHON}" - "${REFERENCE}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree, write_json_exclusive
from scripts.gse105127_exact_reference_contract import admits_exact_reference

reference = Path(sys.argv[1])
output = Path(sys.argv[2])
manifest = verify_frozen_tree(reference)
receipt = json.loads((reference / "receipt.json").read_text(encoding="utf-8"))
if not admits_exact_reference(manifest, receipt):
    raise SystemExit("actual frozen exact-reference contract was not admitted")
write_json_exclusive(
    output / "actual_reference_contract_receipt.json",
    {
        "status": "passed_exact_reference_contract",
        "artifact_class": manifest["metadata"]["artifact_class"],
        "manifest_status": manifest["metadata"]["status"],
        "receipt_status": receipt["status"],
        "labels_accessed": receipt["labels_accessed"],
        "fit_or_score_performed": receipt["fit_or_score_performed"],
    },
)
PY
"${PYTHON}" - <<'PY'
from pathlib import Path

for relative in (
    "scripts/gse105127_exact_reference_contract.py",
    "scripts/build_gse105127_cpg_crosswalk_status_repair.py",
    "scripts/finalize_gse105127_production_activation_status_repair.py",
    "tests/unit/test_gse105127_reference_status_repair.py",
):
    path = Path(relative)
    compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
bash -n slurm/validate_gse105127_reference_status_repair_cpu.sbatch
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/gse105127_exact_reference_contract.py \
  scripts/build_gse105127_cpg_crosswalk_status_repair.py \
  scripts/finalize_gse105127_production_activation_status_repair.py \
  scripts/build_gse105127_cpg_crosswalk.py \
  scripts/finalize_gse105127_production_activation.py \
  tests/unit/test_gse105127_reference_status_repair.py \
  slurm/validate_gse105127_reference_status_repair_cpu.sbatch \
  "${REFERENCE}/ARTIFACTS.json" \
  "${REFERENCE}/receipt.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1])
tests = (root / "unittest.stderr.txt").read_text(encoding="utf-8")
actual = json.loads(
    (root / "actual_reference_contract_receipt.json").read_text(encoding="utf-8")
)
if (
    "Ran 5 tests" not in tests
    or "OK" not in tests
    or "Warning" in tests
    or actual.get("status") != "passed_exact_reference_contract"
    or actual.get("artifact_class") != "gse105127_reference_bundle"
    or actual.get("manifest_status") != "passed_exact_legacy_stream"
    or actual.get("receipt_status") != "passed_exact_legacy_stream"
    or actual.get("labels_accessed") is not False
    or actual.get("fit_or_score_performed") is not False
):
    raise SystemExit("GSE105127 reference-status repair validation differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse105127_reference_status_repair_validation",
        "status": "passed_synthetic_and_actual_contract",
        "locked_v2_sources_modified": False,
        "exact_reference_status_admitted": "passed_exact_legacy_stream",
        "generic_passed_status_admitted": False,
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
