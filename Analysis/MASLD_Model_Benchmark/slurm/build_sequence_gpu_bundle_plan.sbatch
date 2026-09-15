#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-bundle-plan
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
BACKLOG="${ROOT}/config/campaigns/sequence_gpu_dispatcher_backlog_20260824.tsv"
POLICY="${ROOT}/config/campaigns/sequence_gpu_dispatcher_policy_20260824.json"
BUNDLE_001="${ROOT}/config/campaigns/sequence_gpu_bundle_001_20260824.tsv"
BUILDER="${ROOT}/scripts/build_sequence_gpu_bundle_plan.py"
TARGET="${ROOT}/executions/sequence-gpu-bundle-plan-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' 8e858dafe8d0a5055a838fb5e0aadfac9ffe1b8e83e529acec958d29455a2dc6 "${BACKLOG}" | sha256sum --check --strict
printf '%s  %s\n' 36e830550fb11a8ad0293eeab929383968d8dbb66cab67f3ef27a380a899a3de "${POLICY}" | sha256sum --check --strict
printf '%s  %s\n' 042faacf41d698e4003e5ed18829530aa7e9c69b8d3cc5c375775027b7bfcad1 "${BUNDLE_001}" | sha256sum --check --strict
printf '%s  %s\n' 21270672cc097361179eb3f629ee0e607f9872d4edbd17df59b5c5e1b2fb3aaa "${BUILDER}" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
python "${BUILDER}" \
  --backlog "${BACKLOG}" \
  --policy "${POLICY}" \
  --bundle-001 "${BUNDLE_001}" \
  --output "${STAGE}/plan"
python - "${STAGE}/plan" > "${STAGE}/plan/validation.json" <<'PY'
from pathlib import Path
from collections import Counter
import csv
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
with (root / "bundle_allocations.tsv").open(encoding="utf-8", newline="") as handle:
    bundles = list(csv.DictReader(handle, delimiter="\t"))
with (root / "bundle_tasks.tsv").open(encoding="utf-8", newline="") as handle:
    tasks = list(csv.DictReader(handle, delimiter="\t"))
with (root / "excluded_active_or_completed_tasks.tsv").open(encoding="utf-8", newline="") as handle:
    excluded = list(csv.DictReader(handle, delimiter="\t"))
summary = json.loads((root / "bundle_plan_summary.json").read_text(encoding="utf-8"))
if len(bundles) != 33 or len(tasks) != 240 or len(excluded) != 10:
    raise SystemExit("bundled backlog cardinality differs")
if bundles[0]["bundle_id"] != "sequence_gpu_bundle_001":
    raise SystemExit("first durable bundle differs")
if any(not 4 <= int(row["logical_tasks"]) <= 8 for row in bundles):
    raise SystemExit("bundle task count differs")
if any(int(row["requested_gpu_hours"]) > 48 for row in bundles):
    raise SystemExit("bundle wall-time exceeds 48 hours")
if any(row["qos"] != "nslab" or row["partition"] != "gpu" for row in bundles):
    raise SystemExit("future GPU QOS/partition differs")
if any("innovation" in "\t".join(row.values()).lower() for row in bundles + tasks):
    raise SystemExit("innovation routing leaked into future bundle plan")
if any(row["submission_authority"] != "dispatcher_only_no_manual_sbatch" for row in bundles + tasks):
    raise SystemExit("dispatcher-only authority differs")
if Counter(row["logical_work_id"] for row in tasks).most_common(1)[0][1] != 1:
    raise SystemExit("logical task appears in multiple bundles")
if summary["estimated_gpu_hours"] != 930 or summary["requested_gpu_hours"] != 1163:
    raise SystemExit("aggregate resource estimate differs")
print(json.dumps({
    "schema_version": "masld-bench-sequence-gpu-bundle-plan-validation-v1",
    "status": "pass",
    "bundle_allocations": len(bundles),
    "future_logical_tasks": len(tasks),
    "excluded_active_or_completed_tasks": len(excluded),
    "bundle_size_counts": dict(sorted(Counter(row["logical_tasks"] for row in bundles).items())),
    "initial_readiness_counts": dict(sorted(Counter(row["initial_readiness"] for row in bundles).items())),
    "estimated_gpu_hours": summary["estimated_gpu_hours"],
    "requested_gpu_hours": summary["requested_gpu_hours"],
    "gpu_jobs_submitted": 0,
    "future_qos": "nslab",
}, sort_keys=True, indent=2))
PY
sha256sum "${BACKLOG}" "${POLICY}" "${BUNDLE_001}" "${BUILDER}" \
  "${ROOT}/slurm/build_sequence_gpu_bundle_plan.sbatch" > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "plan/bundle_plan_summary.json").read_text())
freeze_tree(root, {
    "artifact_class": "sequence_gpu_bundle_plan",
    "status": summary["status"],
    "bundle_allocations": summary["bundle_allocations"],
    "future_logical_tasks": summary["future_logical_tasks"],
    "excluded_active_or_completed_tasks": summary["excluded_active_or_completed_tasks"],
    "estimated_gpu_hours": summary["estimated_gpu_hours"],
    "requested_gpu_hours": summary["requested_gpu_hours"],
    "manual_gpu_sbatch_prohibited": True,
    "future_qos": "nslab",
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
