#!/usr/bin/env bash
#SBATCH --job-name=model-check-303
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg
export WANDB_MODE=offline
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome/bin/python
CONFIG="${ROOT}/config/alphagenome_gse281364_task_native_activation.json"
SOURCE_ARTIFACT="${ROOT}/executions/model-check-301-21099175"
PRIOR_FAILED="${ROOT}/executions/model-check-302-21099187.failed"
TARGET="${ROOT}/executions/model-check-303-${SLURM_JOB_ID}"
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
  printf 'Refusing to overwrite corrected task-native activation: %s\n' "${TARGET}" >&2
  exit 1
fi

if [[ -d "${PRIOR_FAILED}" && ! -e "${PRIOR_FAILED}/ARTIFACTS.json" ]]; then
  cp --no-clobber "${ROOT}/executions/model-check-302-21099187.out" "${PRIOR_FAILED}/slurm.stdout.txt"
  cp --no-clobber "${ROOT}/executions/model-check-302-21099187.err" "${PRIOR_FAILED}/slurm.stderr.txt"
  "${PYTHON}" - "${PRIOR_FAILED}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


failed = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    failed,
    {
        "artifact_class": "failed_restricted_sequence_task_native_activation",
        "job_id": "21099187",
        "terminal_state": "FAILED",
        "exit_code": "22:0",
        "failure_reason": "anonymous_GitHub_API_request_returned_HTTP_403_before_any_source_was_written",
        "api_prediction_called": False,
        "checkpoint_bytes_downloaded": False,
        "model_forward_executed": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
    },
)
PY
fi

"${PYTHON}" - "${SOURCE_ARTIFACT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree


source = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(source)
PY

install -d -m 0750 "${STAGE}/sources" "${STAGE}/validation"
cp -a "${SOURCE_ARTIFACT}/sources/." "${STAGE}/sources/"

"${PYTHON}" -m unittest -v \
  tests.unit.test_activate_alphagenome_gse281364_task_native \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/activate_alphagenome_gse281364_task_native.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --source-dir "${STAGE}/sources" \
  --output "${STAGE}/activation" \
  > "${STAGE}/validation/activation.stdout.json" \
  2> "${STAGE}/validation/activation.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/activate_alphagenome_gse281364_task_native.py" \
  "${ROOT}/tests/unit/test_activate_alphagenome_gse281364_task_native.py" \
  "${ROOT}/slurm/activate_alphagenome_gse281364_task_native_cpu_v2.sbatch" \
  "${SOURCE_ARTIFACT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
receipt = json.loads((stage / "activation/receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_outcome_blind_task_native_activation_execution_blocked_no_available_credential"
    or receipt["next_allowed_action"]
    != "provide_AlphaGenome_API_key_via_ALPHA_GENOME_API_KEY_without_logging_its_value"
    or receipt["logical_model_version"] != "ALL_FOLDS"
    or receipt["request_manifest"]["requests"] != 1033
    or receipt["request_manifest"]["input_length_bp"] != 524288
    or not receipt["reference_crosswalk"]["all_used_primary_sequence_accessions_equal"]
    or receipt["api_metadata_called"]
    or receipt["api_prediction_called"]
    or receipt["checkpoint_bytes_downloaded"]
    or receipt["model_forward_executed"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_assets_loaded"]
    or receipt["metrics_calculated"]
    or receipt["open_champion_eligible"]
):
    raise SystemExit("corrected AlphaGenome task-native activation receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_sequence_task_native_activation",
        "model_id": "alphagenome",
        "dataset_id": "gse281364",
        "logical_model_version": "ALL_FOLDS",
        "requests": 1033,
        "input_length_bp": 524288,
        "source_snapshot_artifacts_sha256": "50daf76577a08ea811f5a51b37a7e2e0ffea76c92a35fbe95569576ccfd38365",
        "supersedes_for_next_action_string_only": "model-check-301-21099175",
        "api_credential_present": False,
        "huggingface_credential_present": False,
        "api_prediction_called": False,
        "checkpoint_bytes_downloaded": False,
        "model_forward_executed": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "open_champion_eligible": False,
        "status": receipt["status"],
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'Corrected task-native activation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
