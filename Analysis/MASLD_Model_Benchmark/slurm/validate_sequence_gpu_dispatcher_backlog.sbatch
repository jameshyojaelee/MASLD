#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-dispatch-validate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ATTEMPTS="${ROOT}/config/campaigns/sequence_gpu_attempt_registry_20260824.tsv"
BACKLOG="${ROOT}/config/campaigns/sequence_gpu_dispatcher_backlog_20260824.tsv"
POLICY="${ROOT}/config/campaigns/sequence_gpu_dispatcher_policy_20260824.json"
TARGET="${ROOT}/executions/sequence-gpu-dispatcher-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' ecd37495a2ccb0d0077cb5f073683c6765821254259baab06d61f9cfd45357b0 "${ATTEMPTS}" | sha256sum --check --strict
printf '%s  %s\n' 8e858dafe8d0a5055a838fb5e0aadfac9ffe1b8e83e529acec958d29455a2dc6 "${BACKLOG}" | sha256sum --check --strict
printf '%s  %s\n' 36e830550fb11a8ad0293eeab929383968d8dbb66cab67f3ef27a380a899a3de "${POLICY}" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ATTEMPTS}" "${BACKLOG}" "${POLICY}" \
  > "${STAGE}/validation/manifest_validation.json" <<'PY'
from pathlib import Path
from collections import Counter
import csv
import json
import sys

attempt_path, backlog_path, policy_path = map(Path, sys.argv[1:])
with attempt_path.open(encoding="utf-8", newline="") as handle:
    attempts = list(csv.DictReader(handle, delimiter="\t"))
with backlog_path.open(encoding="utf-8", newline="") as handle:
    backlog = list(csv.DictReader(handle, delimiter="\t"))
policy = json.loads(policy_path.read_text(encoding="utf-8"))
logical_ids = {row["logical_work_id"] for row in backlog}
attempt_jobs = {row["slurm_job_id"] for row in attempts}
if len(backlog) != 250 or len(logical_ids) != 250:
    raise SystemExit("logical dispatcher cardinality differs")
if len(attempts) != 232 or len(attempt_jobs) != 232:
    raise SystemExit("attempt registry cardinality differs")
if any(row["logical_work_id"] not in logical_ids for row in attempts):
    raise SystemExit("attempt references unknown logical work")
if Counter(row["model_id"] for row in backlog) != {
    "bpnet": 55,
    "chrombpnet": 55,
    "scbasset": 30,
    "sequence_cnn_control": 55,
    "sequence_transformer_control": 55,
}:
    raise SystemExit("model-family cardinality differs")
if Counter(row["current_state"] for row in backlog) != {
    "active_pending": 1,
    "active_running": 4,
    "backlog_ready": 133,
    "backlog_waiting_input": 92,
    "backlog_waiting_model": 15,
    "completed": 5,
}:
    raise SystemExit("dispatcher state snapshot differs")
normalized_attempt_states = Counter(
    row["slurm_state"].split(" by ", 1)[0] for row in attempts
)
if normalized_attempt_states != {
    "CANCELLED": 218,
    "COMPLETED": 5,
    "FAILED": 4,
    "PENDING": 1,
    "RUNNING": 4,
}:
    raise SystemExit("attempt state snapshot differs")
active = {
    row["active_job_id"]
    for row in backlog
    if row["current_state"] in {"active_running", "active_pending"}
}
if active != {"21066983", "21066987", "21066988", "21066995", "21067819"}:
    raise SystemExit("active dispatcher identity differs")
if any(
    row["submission_authority"] != "dispatcher_only_no_manual_sbatch"
    or row["routing_policy"] != "nslab_only"
    or row["gpu_squeue_ceiling"] != "5"
    or row["max_running_gpu"] != "4"
    or row["max_pending_gpu"] != "1"
    or "outcome" in row["export_contract"].lower()
    or "sealed" in row["export_contract"].lower()
    for row in backlog
):
    raise SystemExit("dispatcher safety contract differs")
for row in backlog:
    prerequisite = row["prerequisite_logical_id"]
    if row["stage"] == "predict_valid":
        if prerequisite not in logical_ids or "__train__" not in prerequisite:
            raise SystemExit("scBasset prediction dependency differs")
    elif prerequisite != "none":
        raise SystemExit("unexpected model prerequisite")
if (
    policy["gpu_squeue_ceiling"] != 5
    or policy["max_running_gpu_jobs"] != 4
    or policy["max_pending_gpu_jobs"] != 1
    or policy["submission_batch_size"] != 1
    or not policy["manual_gpu_sbatch_prohibited"]
    or policy["routing"] != {
        "account": "nslab",
        "qos": "nslab",
        "gres": "gpu:l40s:1",
        "policy": "nslab_only",
    }
):
    raise SystemExit("dispatcher policy differs")
print(json.dumps({
    "schema_version": "masld-bench-sequence-gpu-dispatcher-validation-v1",
    "status": "pass",
    "logical_tasks": len(backlog),
    "attempts": len(attempts),
    "active_gpu_jobs": sorted(active),
    "state_counts": dict(sorted(Counter(row["current_state"] for row in backlog).items())),
    "attempt_state_counts": dict(sorted(normalized_attempt_states.items())),
    "manual_gpu_submission_prohibited": True,
}, sort_keys=True, indent=2))
PY
squeue -u "${USER}" -h -o '%i|%t|%b' \
  | awk -F '|' '$3 ~ /gpu/ {print}' > "${STAGE}/validation/live_gpu_squeue.tsv"
python - "${STAGE}/validation/live_gpu_squeue.tsv" \
  > "${STAGE}/validation/live_gpu_guard.json" <<'PY'
from pathlib import Path
from collections import Counter
import json
import sys

rows = [line.split("|") for line in Path(sys.argv[1]).read_text().splitlines() if line]
states = Counter(row[1] for row in rows)
if len(rows) > 5 or states["R"] > 4 or states["PD"] > 1:
    raise SystemExit("live GPU queue exceeds dispatcher ceiling")
print(json.dumps({
    "status": "pass", "gpu_jobs": len(rows),
    "running": states["R"], "pending": states["PD"],
}, sort_keys=True, indent=2))
PY
sha256sum "${ATTEMPTS}" "${BACKLOG}" "${POLICY}" \
  "${ROOT}/slurm/validate_sequence_gpu_dispatcher_backlog.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
manifest = json.loads((root / "validation/manifest_validation.json").read_text())
guard = json.loads((root / "validation/live_gpu_guard.json").read_text())
freeze_tree(root, {
    "artifact_class": "sequence_gpu_dispatcher_validation",
    "status": "passed",
    "logical_tasks": manifest["logical_tasks"],
    "attempts": manifest["attempts"],
    "gpu_jobs_at_validation": guard["gpu_jobs"],
    "running_gpu_jobs_at_validation": guard["running"],
    "pending_gpu_jobs_at_validation": guard["pending"],
    "manual_gpu_submission_prohibited": True,
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
