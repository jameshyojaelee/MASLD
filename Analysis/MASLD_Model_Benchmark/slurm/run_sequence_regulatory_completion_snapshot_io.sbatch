#!/usr/bin/env bash
#SBATCH --job-name=model-audit
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-audit-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-audit-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT=${ROOT}/config/evaluation/sequence_regulatory_full_rectangle_completion_20260825.json
AUDIT=${ROOT}/scripts/audit_sequence_regulatory_full_rectangle_completion.py
ACCEPTANCE=${ROOT}/executions/sequence-regulatory-completion-gate-tests-21100975
TARGET=${ROOT}/executions/sequence-regulatory-completion-snapshot-${SLURM_JOB_ID}
STAGE=${TARGET}.staging

preserve_failed_stage() {
  local failed=${TARGET}.failed
  if [[ -d "${STAGE}" && ! -e "${failed}" ]]; then
    mv "${STAGE}" "${failed}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite completion snapshot: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' \
  edede235b0b87913de334142a25d184a89803870cb92c640a54a81375d68ec08 "${CONTRACT}" \
  afba2e8882c84596451a0d74ac387b6ca587bdbde825c35352b5246ade8ca491 "${AUDIT}" \
  4632e8d6debea40b832365c31fffd301bcd9896dd765f2722f3b2332ff1e4d29 "${ACCEPTANCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
PYTHONPATH=${ROOT}/src:${ROOT} python - "${ACCEPTANCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
PY

PYTHONPATH=${ROOT}/src:${ROOT} python "${AUDIT}" \
  --root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/completion" \
  > "${STAGE}/validation/audit.stdout.json"
python --version > "${STAGE}/validation/python-version.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum "${CONTRACT}" "${AUDIT}" "$0" > "${STAGE}/source.sha256"

PYTHONPATH=${ROOT}/src:${ROOT} python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads(
    (root / "completion/completion_lock.json").read_text(encoding="utf-8")
)
if (
    receipt["status"] not in {
        "waiting_for_full_rectangles",
        "blocked_rectangle_requires_new_immutable_retry",
        "ready_both_rectangles_locked",
    }
    or receipt["partial_ranking_authorized"]
    or receipt["cross_task_family_ranking_authorized"]
    or receipt["prediction_values_read"]
    or receipt["metric_values_read"]
    or receipt["observed_outcomes_read"]
    or receipt["sealed_data_read"]
    or receipt["champion_claim_allowed"]
    or receipt["universal_claim_allowed"]
):
    raise SystemExit("completion snapshot firewall differs")
print(freeze_tree(root, {
    "artifact_class": "sequence_regulatory_completion_snapshot",
    "status": receipt["status"],
    "both_full_rectangles": receipt["both_full_rectangles"],
    "evaluation_launch_authorized": receipt["evaluation_launch_authorized"],
    "sequence_complete_fits": receipt["task_native_sequence"]["fit_disposition_counts"].get("complete_compatible", 0),
    "sequence_complete_prediction_views": receipt["task_native_sequence"]["prediction_view_disposition_counts"].get("complete_compatible", 0),
    "scbasset_fit_disposition_counts": receipt["scbasset"]["fit_disposition_counts"],
    "outcomes_read": False,
    "prediction_values_read": False,
    "sealed_data_read": False,
    "partial_ranking_performed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Completion snapshot: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
