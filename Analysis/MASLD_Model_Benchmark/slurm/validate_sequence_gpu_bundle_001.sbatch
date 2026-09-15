#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-bundle001-check
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
MANIFEST="${ROOT}/config/campaigns/sequence_gpu_bundle_001_20260824.tsv"
BACKLOG="${ROOT}/config/campaigns/sequence_gpu_dispatcher_backlog_20260824.tsv"
WRAPPER="${ROOT}/slurm/sequence_gpu_bundle_001.sbatch"
TARGET="${ROOT}/executions/sequence-gpu-bundle-001-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 042faacf41d698e4003e5ed18829530aa7e9c69b8d3cc5c375775027b7bfcad1 "${MANIFEST}" | sha256sum --check --strict
printf '%s  %s\n' 8e858dafe8d0a5055a838fb5e0aadfac9ffe1b8e83e529acec958d29455a2dc6 "${BACKLOG}" | sha256sum --check --strict
printf '%s  %s\n' c2ca6a2d2cfce730c55a767f4cadf4b5d846de2be90d8b51380852f607fd3d05 "${WRAPPER}" | sha256sum --check --strict
bash -n "${WRAPPER}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${MANIFEST}" "${BACKLOG}" "${WRAPPER}" \
  > "${STAGE}/validation/bundle_validation.json" <<'PY'
from pathlib import Path
from collections import Counter
import csv
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

manifest_path, backlog_path, wrapper_path = map(Path, sys.argv[1:])
with manifest_path.open(encoding="utf-8", newline="") as handle:
    tasks = list(csv.DictReader(handle, delimiter="\t"))
with backlog_path.open(encoding="utf-8", newline="") as handle:
    backlog = {row["logical_work_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
if len(tasks) != 4 or len({row["logical_work_id"] for row in tasks}) != 4:
    raise SystemExit("bundle 001 cardinality differs")
if Counter(row["model_id"] for row in tasks) != {"bpnet": 2, "chrombpnet": 2}:
    raise SystemExit("bundle 001 model balance differs")
if Counter(row["lineage_id"] for row in tasks) != {"cholangiocyte": 2, "fibroblast": 2}:
    raise SystemExit("bundle 001 lineage balance differs")
for order, task in enumerate(tasks, 1):
    if int(task["task_order"]) != order:
        raise SystemExit("bundle 001 task order differs")
    original = backlog.get(task["logical_work_id"])
    if original is None or original["current_state"] != "backlog_ready":
        raise SystemExit("bundle task is not backlog-ready")
    for key in ("model_id", "stage", "lineage_id", "outer_fold", "seed",
                "input_job_id", "input_artifacts_sha256", "script"):
        if task[key] != original[key]:
            raise SystemExit(f"bundle binding differs for {key}")
    if task["prerequisite_logical_id"] != "none" or task["failure_scope"] != "isolated_continue":
        raise SystemExit("bundle failure/dependency policy differs")
    root = manifest_path.parents[2]
    split_id = f"donor{task['outer_fold']}_genomic{task['outer_fold']}"
    artifact = root / "executions" / f"chrombpnet-{task['lineage_id']}-{split_id}-inputs-{task['input_job_id']}"
    verify_frozen_tree(artifact)
    import hashlib
    actual = hashlib.sha256((artifact / "ARTIFACTS.json").read_bytes()).hexdigest()
    if actual != task["input_artifacts_sha256"]:
        raise SystemExit("bundle input artifact hash differs")
wrapper = wrapper_path.read_text(encoding="utf-8")
required = [
    "#SBATCH --gres=gpu:l40s:1",
    "#SBATCH --cpus-per-task=8",
    "#SBATCH --mem=32G",
    "#SBATCH --time=20:00:00",
    "MASLD_GPU_DISPATCHER_AUTHORITY",
    "isolated_continue",
    "task_receipts",
]
if any(value not in wrapper for value in required):
    raise SystemExit("bundle wrapper contract differs")
print(json.dumps({
    "schema_version": "masld-bench-sequence-gpu-bundle-validation-v1",
    "status": "pass",
    "bundle_id": "sequence_gpu_bundle_001",
    "logical_tasks": len(tasks),
    "models": dict(sorted(Counter(row["model_id"] for row in tasks).items())),
    "estimated_gpu_hours": 16,
    "requested_gpu_hours": 20,
    "gpu_jobs_submitted": 0,
    "dispatcher_authority_required": True,
}, sort_keys=True, indent=2))
PY
sha256sum "${MANIFEST}" "${BACKLOG}" "${WRAPPER}" \
  "${ROOT}/slurm/validate_sequence_gpu_bundle_001.sbatch" > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "validation/bundle_validation.json").read_text())
freeze_tree(root, {
    "artifact_class": "sequence_gpu_bundle_validation",
    "bundle_id": summary["bundle_id"],
    "status": summary["status"],
    "logical_tasks": summary["logical_tasks"],
    "estimated_gpu_hours": summary["estimated_gpu_hours"],
    "requested_gpu_hours": summary["requested_gpu_hours"],
    "gpu_jobs_submitted": 0,
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
