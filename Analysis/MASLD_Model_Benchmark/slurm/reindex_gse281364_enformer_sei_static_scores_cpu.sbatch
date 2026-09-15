#!/usr/bin/env bash
#SBATCH --job-name=model-check-274
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_enformer_sei_static_reindex.json
TARGET=${ROOT}/executions/model-check-274-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Enformer/Sei reindex attempt: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src
cd "${ROOT}"

"${PYTHON}" -m unittest -v \
  tests.unit.test_reindex_gse281364_enformer_sei_static_scores \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/reindex_gse281364_enformer_sei_static_scores.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/binding" \
  > "${STAGE}/validation/reindex.stdout.json" \
  2> "${STAGE}/validation/reindex.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/reindex_gse281364_enformer_sei_static_scores.py" \
  "${ROOT}/tests/unit/test_reindex_gse281364_enformer_sei_static_scores.py" \
  "${ROOT}/slurm/reindex_gse281364_enformer_sei_static_scores_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "binding/receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_outcome_blind_static_reindex"
    or receipt["elements"] != 1033
    or receipt["row_bindings"] != 10330
    or receipt["long_range_blocks"] != 239
    or receipt["source_groups_reassigned_from_legacy_folds"] != 830
    or receipt["models"] != ["enformer_crested_restricted_port", "sei"]
    or not receipt["prediction_values_read"]
    or receipt["outcomes_read"]
    or receipt["sealed_assets_read"]
    or receipt["model_fit"]
    or receipt["calibration_fit"]
    or receipt["models_ranked"]
    or receipt["native_outputs_are_interchangeable"]
    or receipt["global_frozen_census_changed"]
    or receipt["open_champion_eligible"]
):
    raise SystemExit("Enformer/Sei static reindex receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "gse281364_enformer_sei_static_score_reindex",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "models": ["enformer_crested_restricted_port", "sei"],
        "elements": 1033,
        "row_bindings": 10330,
        "long_range_blocks": 239,
        "prediction_values_read": True,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "model_fit": False,
        "models_ranked": False,
        "global_frozen_census_changed": False,
        "status": "passed_outcome_blind_static_reindex",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Enformer/Sei static-score reindex: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
