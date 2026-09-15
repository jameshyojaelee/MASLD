#!/usr/bin/env bash
#SBATCH --job-name=model-check-ep904
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
ADMISSION="${ROOT}/executions/epibert-epcotv2-admission-sources-21064417"
RUNTIME="${ROOT}/executions/environments/epibert-executable-tf2.13-21064522"
PRIMARY="${ROOT}/executions/model-data-ep903-21098560"
TARGET="${ROOT}/executions/model-check-ep904-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
export PYTHONPATH="${ROOT}:${ROOT}/src"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite EpiBERT native preprocessing audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_epibert_native_preprocessing \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/audit_epibert_native_preprocessing.py" \
  --admission "${ADMISSION}" \
  --runtime "${RUNTIME}" \
  --primary-evidence "${PRIMARY}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
/usr/bin/git --version > "${STAGE}/validation/git.version.txt"
sha256sum \
  "${ROOT}/scripts/audit_epibert_native_preprocessing.py" \
  "${ROOT}/tests/unit/test_audit_epibert_native_preprocessing.py" \
  "${ROOT}/slurm/audit_epibert_native_preprocessing_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
receipt = json.loads(
    (stage / "audit/native_preprocessing_audit.json").read_text(encoding="utf-8")
)
tests = (stage / "validation/tests.stderr.txt").read_text(encoding="utf-8")
required_unresolved = {
    "atac_scale_factor_units",
    "atac_scale_numerator",
    "motif_top_peaks",
    "motif_peak_window",
    "motif_background_peaks",
    "motif_sea_threshold",
    "motif_normalization",
    "peak_calling_and_reference",
}
if (
    receipt["status"] != "pass_terminal_blocker_frozen"
    or receipt["native_preprocessing_exact"]
    or receipt["native_preprocessing_fixture_built"]
    or receipt["terminal_disposition"]
    != "checkpoint_runtime_ready_native_preprocessing_terminally_blocked_no_fixture"
    or receipt["field_disposition"]["motif_roster_and_order"]["count"] != 693
    or receipt["field_disposition"]["motif_roster_and_order"]["status"] != "authoritative"
    or not required_unresolved.issubset(set(receipt["unresolved_fields"]))
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["sealed_assets_read"]
    or receipt["metrics_calculated"]
    or receipt["gpu_queue_item_created"]
    or "Ran 2 tests" not in tests
    or "OK" not in tests
):
    raise SystemExit("EpiBERT native preprocessing audit differs")
freeze_tree(
    stage,
    {
        "artifact_class": "epibert_native_preprocessing_terminal_audit",
        "model_id": "epibert",
        "checkpoint_revision": "dae61b434f885991c29e2fde3fbf7f70d0ff9b80",
        "primary_evidence_artifacts_sha256": "4047b17e387eaf1be1401a7e80bb9470670ef9af11b13e2f4209a25133e1e6e0",
        "native_preprocessing_exact": False,
        "fixture_built": False,
        "task_scope": "observed_ATAC_only",
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "status": "pass_terminal_blocker_frozen",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'EpiBERT native preprocessing audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
