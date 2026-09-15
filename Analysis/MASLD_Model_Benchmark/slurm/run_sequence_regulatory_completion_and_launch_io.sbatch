#!/usr/bin/env bash
#SBATCH --job-name=model-audit
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
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
WRAPPER_ACCEPTANCE=${ROOT}/executions/sequence-regulatory-postprocess-wrapper-tests-21101011
SEQUENCE_EVALUATOR=${ROOT}/slurm/evaluate_sequence_task_native_five_seed_rectangle_cpu.sbatch
SCBASSET_EVALUATOR=${ROOT}/slurm/evaluate_scbasset_five_seed_rectangle_cpu.sbatch
TARGET=${ROOT}/executions/sequence-regulatory-completion-lock-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
LAUNCH_TARGET=${ROOT}/executions/sequence-regulatory-evaluation-launch-${SLURM_JOB_ID}
LAUNCH_STAGE=${LAUNCH_TARGET}.staging

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${TARGET}.failed" ]]; then
    mv "${STAGE}" "${TARGET}.failed"
  fi
  if [[ -d "${LAUNCH_STAGE}" && ! -e "${LAUNCH_TARGET}.failed" ]]; then
    mv "${LAUNCH_STAGE}" "${LAUNCH_TARGET}.failed"
  fi
}
trap preserve_failed_stage EXIT
for value in "${TARGET}" "${STAGE}" "${TARGET}.failed" \
  "${LAUNCH_TARGET}" "${LAUNCH_STAGE}" "${LAUNCH_TARGET}.failed"; do
  if [[ -e "${value}" ]]; then
    printf 'Refusing to overwrite production gate artifact: %s\n' "${value}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' \
  edede235b0b87913de334142a25d184a89803870cb92c640a54a81375d68ec08 "${CONTRACT}" \
  afba2e8882c84596451a0d74ac387b6ca587bdbde825c35352b5246ade8ca491 "${AUDIT}" \
  4632e8d6debea40b832365c31fffd301bcd9896dd765f2722f3b2332ff1e4d29 "${ACCEPTANCE}/ARTIFACTS.json" \
  2359c0bb1e075809787c6162f6539c61b1e300d1659f44996dc291e77755f676 "${WRAPPER_ACCEPTANCE}/ARTIFACTS.json" \
  c134bd7796af7b28d75e4a9cfbe4f75a045413d14e86695a6309e538f6482611 "${SEQUENCE_EVALUATOR}" \
  41d1c755497d105b979f7cca7349a0af936ef11e1ae5488dd3127d73f186392c "${SCBASSET_EVALUATOR}" \
  | sha256sum --check --strict
PYTHONPATH=${ROOT}/src:${ROOT} python - "${ACCEPTANCE}" "${WRAPPER_ACCEPTANCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY

PYTHONPATH=${ROOT}/src:${ROOT} python "${AUDIT}" \
  --root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/completion" \
  > "${STAGE}/validation/audit.stdout.json"
python --version > "${STAGE}/validation/python-version.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONTRACT}" "${AUDIT}" "${SEQUENCE_EVALUATOR}" \
  "${SCBASSET_EVALUATOR}" "$0" > "${STAGE}/source.sha256"

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
    raise SystemExit("production completion firewall differs")
print(freeze_tree(root, {
    "artifact_class": "sequence_regulatory_completion_lock",
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
LOCK_SHA=$(sha256sum "${TARGET}/ARTIFACTS.json" | awk '{print $1}')
STATUS=$(python - "${TARGET}/completion/completion_lock.json" <<'PY'
from pathlib import Path
import json
import sys
print(json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))["status"])
PY
)
if [[ "${STATUS}" != ready_both_rectangles_locked ]]; then
  trap - EXIT
  printf 'Completion gate terminal without evaluation launch: %s %s\n' \
    "${STATUS}" "${TARGET}"
  printf '%s  %s\n' "${LOCK_SHA}" "${TARGET}/ARTIFACTS.json"
  exit 0
fi

PRIOR=$(PYTHONPATH=${ROOT}/src:${ROOT} python - "${ROOT}" "${LOCK_SHA}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
target = sys.argv[2]
matches = []
for path in sorted((root / "executions").glob("sequence-regulatory-evaluation-launch-*/launch.json")):
    try:
        row = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        continue
    if row.get("completion_artifacts_sha256") == target:
        matches.append(str(path.parent))
print("\n".join(matches))
PY
)
if [[ -n "${PRIOR}" ]]; then
  printf 'Evaluation already launched for completion %s: %s\n' "${LOCK_SHA}" "${PRIOR}" >&2
  exit 1
fi

SEQUENCE_JOB=$(sbatch --parsable \
  --export=ALL,COMPLETION_ROOT="${TARGET}",COMPLETION_SHA="${LOCK_SHA}" \
  "${SEQUENCE_EVALUATOR}")
SCBASSET_JOB=$(sbatch --parsable \
  --export=ALL,COMPLETION_ROOT="${TARGET}",COMPLETION_SHA="${LOCK_SHA}" \
  "${SCBASSET_EVALUATOR}")
install -d -m 0750 "${LAUNCH_STAGE}/validation"
python - "${LAUNCH_STAGE}/launch.json" "${TARGET}" "${LOCK_SHA}" \
  "${SEQUENCE_JOB}" "${SCBASSET_JOB}" <<'PY'
from pathlib import Path
import json
import sys

output, completion, digest, sequence_job, scbasset_job = sys.argv[1:]
Path(output).write_text(json.dumps({
    "schema_version": "masld-bench-sequence-regulatory-evaluation-launch-v1",
    "status": "submitted",
    "completion_root": completion,
    "completion_artifacts_sha256": digest,
    "task_native_sequence_evaluation_job_id": int(sequence_job),
    "scbasset_evaluation_job_id": int(scbasset_job),
    "partial_ranking_performed": False,
    "cross_task_family_ranking_performed": False,
    "sealed_data_read": False,
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
sha256sum "${SEQUENCE_EVALUATOR}" "${SCBASSET_EVALUATOR}" "$0" \
  > "${LAUNCH_STAGE}/source.sha256"
PYTHONPATH=${ROOT}/src:${ROOT} python - "${LAUNCH_STAGE}" "${LOCK_SHA}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(freeze_tree(root, {
    "artifact_class": "sequence_regulatory_evaluation_launch",
    "status": "submitted",
    "completion_artifacts_sha256": sys.argv[2],
    "evaluation_jobs": 2,
    "partial_ranking_performed": False,
    "cross_task_family_ranking_performed": False,
    "sealed_data_read": False,
}))
PY
mv "${LAUNCH_STAGE}" "${LAUNCH_TARGET}"
trap - EXIT
printf 'Completion lock: %s %s\n' "${TARGET}" "${LOCK_SHA}"
printf 'Evaluation jobs: task-native=%s scBasset=%s\n' \
  "${SEQUENCE_JOB}" "${SCBASSET_JOB}"
sha256sum "${LAUNCH_TARGET}/ARTIFACTS.json"
