#!/usr/bin/env bash
#SBATCH --job-name=model-data-928
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-880-21130257-gse135251-source
SOURCE_ARTIFACTS_SHA256=4dab6a759371ccca56227bc3ce1df2c4e9ad27397cab07db4c337669c5003e8d
TARGET=${ROOT}/executions/model-data-064-21079902/fixture/molecular
TARGET_ARTIFACTS_SHA256=2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456
JOIN=${ROOT}/executions/gse267145-authoritative-join-21064930
JOIN_ARTIFACTS_SHA256=e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580
OUT=${ROOT}/executions/model-data-928-${SLURM_JOB_ID}
STAGE=${OUT}.staging
FAILED=${OUT}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${OUT}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite contamination-audit artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${TARGET_ARTIFACTS_SHA256}" "${TARGET}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${JOIN_ARTIFACTS_SHA256}" "${JOIN}/ARTIFACTS.json" | sha256sum --check --strict
bash -n slurm/audit_gse135251_gse267145_contamination_cpu.sbatch
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_gse135251_gse267145_contamination \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/audit_gse135251_gse267145_contamination.py \
  --benchmark-root "${ROOT}" \
  --source-root "${SOURCE}/source" \
  --target-root "${TARGET}" \
  --join-root "${JOIN}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/audit_stdout.json" 2> "${STAGE}/audit_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/audit_gse135251_gse267145_contamination.py \
  tests/unit/test_audit_gse135251_gse267145_contamination.py \
  slurm/audit_gse135251_gse267145_contamination_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${TARGET}/ARTIFACTS.json" \
  "${JOIN}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = root / "audit"
receipt = json.loads((audit / "audit.json").read_text(encoding="utf-8"))
if receipt.get("status") not in {"pass_no_overlap_detected", "fail_overlap_detected"}:
    raise SystemExit("contamination audit status differs")
if receipt["identifier_layer"]["shared_identifiers_total"] != 0:
    raise SystemExit("identifier overlap detected between GSE135251 and GSE267145")
if receipt["molecular_layer"]["reciprocal_best_hits_exceeding_reference"]:
    raise SystemExit("molecular near-duplicate detected between the two cohorts")
if not receipt["identifier_layer"]["namespaces_compared"]:
    raise SystemExit("no identifier namespace was actually comparable")
freeze_tree(audit, {
    "artifact_class": "gse135251_gse267145_contamination_audit",
    "source_series": "GSE135251",
    "target_series": "GSE267145",
    "identifier_overlap": 0,
    "molecular_near_duplicates": 0,
    "status": receipt["status"],
})
verify_frozen_tree(audit)
freeze_tree(root, {
    "artifact_class": "gse135251_gse267145_contamination_audit_campaign",
    "status": receipt["status"],
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${OUT}"
sha256sum "${OUT}/ARTIFACTS.json" "${OUT}/audit/ARTIFACTS.json"
