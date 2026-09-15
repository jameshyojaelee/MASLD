#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-bundle23-check
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PLAN="${ROOT}/executions/sequence-gpu-bundle-plan-21068430"
CONTROLLER="${ROOT}/scripts/run_sequence_gpu_bundle.py"
WRAPPER="${ROOT}/slurm/run_sequence_gpu_bundle.sbatch"
FAILED_ATTEMPTS="${ROOT}/config/campaigns/sequence_gpu_bundle_failed_admission_attempts_20260824.tsv"
TARGET="${ROOT}/executions/sequence-gpu-bundles-002-003-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 \
  "${STAGE}/validation/bundle_002/logs" \
  "${STAGE}/validation/bundle_002/task_receipts" \
  "${STAGE}/validation/bundle_003/logs" \
  "${STAGE}/validation/bundle_003/task_receipts"
printf '%s  %s\n' 36dd9ad437bdd61e150e1c681a88bfdf946672c768bffe5f78fac8ab2c352e69 "${PLAN}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' fa31c4cbb5c7ff49746c19d9a3bb49621cdf3205f9e485c17952e3c567ec078f "${CONTROLLER}" | sha256sum --check --strict
printf '%s  %s\n' b4b27edddfac127a50a5065bc96896adbe6a77048d01a392a8577558dd2b46c0 "${WRAPPER}" | sha256sum --check --strict
printf '%s  %s\n' 06c9622e88bdb32d6387ddb8bc540d22aa86ade241f72f79c419a4016c1ae140 "${FAILED_ATTEMPTS}" | sha256sum --check --strict
bash -n "${WRAPPER}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PLAN}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
for bundle_number in 002 003; do
  env \
    MASLD_GPU_DISPATCHER_AUTHORITY=sequence_gpu_dispatcher_v2 \
    BUNDLE_PLAN_ROOT="${PLAN}" \
    BUNDLE_PLAN_ARTIFACTS_SHA256=36dd9ad437bdd61e150e1c681a88bfdf946672c768bffe5f78fac8ab2c352e69 \
    BUNDLE_ID="sequence_gpu_bundle_${bundle_number}" \
    BUNDLE_REQUESTED_TIME=40:00:00 \
    BUNDLE_PREFLIGHT_ONLY=1 \
    BUNDLE_PREFLIGHT_STAGE="${STAGE}/validation/bundle_${bundle_number}" \
    SLURM_JOB_ID="990000${bundle_number}" \
    bash "${WRAPPER}"
done
python - "${PLAN}" "${STAGE}/validation" "${WRAPPER}" "${FAILED_ATTEMPTS}" \
  > "${STAGE}/validation/summary.json" <<'PY'
from pathlib import Path
from collections import Counter
import csv
import json
import sys

plan, validation, wrapper_path, failed_path = map(Path, sys.argv[1:])
with (plan / "plan/bundle_allocations.tsv").open(encoding="utf-8", newline="") as handle:
    allocations = {row["bundle_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
with (plan / "plan/bundle_tasks.tsv").open(encoding="utf-8", newline="") as handle:
    all_tasks = list(csv.DictReader(handle, delimiter="\t"))
seen = set()
for number in ("002", "003"):
    bundle_id = f"sequence_gpu_bundle_{number}"
    allocation = allocations[bundle_id]
    tasks = [row for row in all_tasks if row["bundle_id"] == bundle_id]
    if (
        len(tasks) != 8
        or allocation["models"] != "bpnet,chrombpnet"
        or allocation["initial_readiness"] != "ready"
        or allocation["requested_time"] != "40:00:00"
        or allocation["qos"] != "nslab"
        or allocation["partition"] != "gpu"
    ):
        raise SystemExit(f"{bundle_id} allocation contract differs")
    if Counter(row["model_id"] for row in tasks) != {"bpnet": 4, "chrombpnet": 4}:
        raise SystemExit(f"{bundle_id} model balance differs")
    if Counter(row["lineage_id"] for row in tasks) != {
        "cholangiocyte": 2, "fibroblast": 2, "macrophage": 2, "t_cell": 2,
    }:
        raise SystemExit(f"{bundle_id} lineage balance differs")
    if any(row["current_state"] != "backlog_ready" for row in tasks):
        raise SystemExit(f"{bundle_id} contains a non-ready task")
    logical_ids = {row["logical_work_id"] for row in tasks}
    if seen & logical_ids:
        raise SystemExit("bundle 002/003 logical tasks overlap")
    seen |= logical_ids
    preflight = json.loads(
        (validation / f"bundle_{number}/preflight.json").read_text(encoding="utf-8")
    )
    if (
        preflight["status"] != "pass"
        or preflight["logical_tasks"] != 8
        or preflight["outcome_paths_available_to_models"]
        or preflight["sealed_paths_available_to_models"]
    ):
        raise SystemExit(f"{bundle_id} preflight differs")
wrapper = wrapper_path.read_text(encoding="utf-8")
for required in (
    "#SBATCH --account=nslab", "#SBATCH --qos=nslab",
    "#SBATCH --gres=gpu:l40s:1", "#SBATCH --time=48:00:00",
    "MASLD_GPU_DISPATCHER_AUTHORITY", "BUNDLE_REQUESTED_TIME",
):
    if required not in wrapper:
        raise SystemExit("generic bundle wrapper contract differs")
if "innovation" in wrapper.lower():
    raise SystemExit("innovation routing present in generic bundle wrapper")
with failed_path.open(encoding="utf-8", newline="") as handle:
    failed = list(csv.DictReader(handle, delimiter="\t"))
if (
    [(row["slurm_job_id"], row["bundle_id"]) for row in failed]
    != [("21068456", "sequence_gpu_bundle_002"), ("21068461", "sequence_gpu_bundle_003")]
    or any(
        row["slurm_state"] != "FAILED"
        or row["failure_stage"] != "wrapper_admission"
        or row["model_code_started"] != "no"
        or row["input_data_read"] != "no"
        for row in failed
    )
):
    raise SystemExit("failed admission-attempt registry differs")
print(json.dumps({
    "schema_version": "masld-bench-sequence-gpu-bundle-validation-v2",
    "status": "pass",
    "bundles": ["sequence_gpu_bundle_002", "sequence_gpu_bundle_003"],
    "logical_tasks": len(seen),
    "model_task_counts": {"bpnet": 8, "chrombpnet": 8},
    "requested_gpu_hours_per_bundle": 40,
    "gpu_jobs_submitted": 0,
    "future_qos": "nslab",
    "dispatcher_authority_required": True,
    "preserved_failed_attempts": ["21068456", "21068461"],
    "wrapper_preflight_invoked_for_exact_plan_ids": True,
}, sort_keys=True, indent=2))
PY
sha256sum "${PLAN}/ARTIFACTS.json" "${CONTROLLER}" "${WRAPPER}" "${FAILED_ATTEMPTS}" \
  "${ROOT}/slurm/validate_sequence_gpu_bundles_002_003.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "validation/summary.json").read_text())
freeze_tree(root, {
    "artifact_class": "sequence_gpu_bundle_validation",
    "bundles": summary["bundles"],
    "status": summary["status"],
    "logical_tasks": summary["logical_tasks"],
    "model_task_counts": summary["model_task_counts"],
    "gpu_jobs_submitted": 0,
    "future_qos": "nslab",
    "preserved_failed_attempts": summary["preserved_failed_attempts"],
    "wrapper_preflight_invoked_for_exact_plan_ids": True,
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
