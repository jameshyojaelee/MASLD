#!/usr/bin/env bash
#SBATCH --job-name=model-data-283
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
CONFIG="${ROOT}/config/scooby_epicardioids_native_context_release.json"
TARGET="${ROOT}/executions/model-data-283-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Scooby native context release: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_acquire_scooby_native_context_release \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" "${ROOT}/scripts/acquire_scooby_native_context_release.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/context" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"

"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/acquire_scooby_native_context_release.py" \
  "${ROOT}/tests/unit/test_acquire_scooby_native_context_release.py" \
  "${ROOT}/slurm/acquire_scooby_native_context_release_io.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
receipt = json.loads((stage / "context/receipt.json").read_text(encoding="utf-8"))
boundary = receipt["release_boundary"]
if (
    receipt["status"] != "pass_native_context_rows_released_query_encoder_still_absent"
    or len(receipt["members"]) != 2
    or receipt["context_encoder_weight_candidates"]
    or not receipt["native_context_rows_downloaded"]
    or receipt["native_context_rows_deserialized"]
    or receipt["full_archive_downloaded"]
    or receipt["checkpoint_loaded"]
    or receipt["model_forward_executed"]
    or receipt["project_data_read"]
    or receipt["sealed_data_read"]
    or boundary["exact_scglue_encoder_weights_released_in_archive"]
    or boundary["coordinate_compatible_liver_query_supported"]
    or boundary["rna_conditioned_atac_eligible"]
):
    raise SystemExit("Scooby native context release receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_observed_atac_model_native_context_release",
        "model_id": "scooby_epicardioids",
        "native_context_rows_downloaded": True,
        "native_context_rows_deserialized": False,
        "exact_scglue_query_encoder_available": False,
        "full_archive_downloaded": False,
        "model_forward_executed": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "rna_conditioned_atac_eligible": False,
        "open_champion_eligible": False,
        "status": "pass_native_context_rows_released_query_encoder_still_absent",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'Scooby native context release: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
