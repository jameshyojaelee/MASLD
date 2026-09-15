#!/usr/bin/env bash
#SBATCH --job-name=model-check-339
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
STAGED="${ROOT}/config/campaigns/gpu_bundle_queue_staging/100-model-training-333.json"
LIVE="${ROOT}/config/campaigns/gpu_bundle_queue/100-model-training-333.json"
QUEUE_SHA256=be194e31dac2e9e08ebef934c79638dfbf118374cd335b997fad05b722bdca89
TARGET="${ROOT}/executions/peakvi-cross-cohort-revision2-queue-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${LIVE}" ]]; then
  printf 'PeakVI live queue target or admission output already exists\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${QUEUE_SHA256}" "${STAGED}" | sha256sum --check --strict
cd "${ROOT}"
PYTHONPATH="src:." "${PYTHON}" -m unittest -v \
  tests.unit.test_gpu_bundle_dispatcher \
  tests.unit.test_gpu_dependency_dispatcher \
  > "${STAGE}/tests.stdout.txt" \
  2> "${STAGE}/tests.stderr.txt"
sha256sum \
  "${STAGED}" \
  config/campaigns/peakvi_gse296875_to_gse281367_fit_predict_revision2_20260825.json \
  scripts/gpu_bundle_dispatcher.py \
  scripts/gpu_dependency_dispatcher.py \
  slurm/run_peakvi_cross_cohort_revision2_bundle.sbatch \
  slurm/validate_peakvi_cross_cohort_revision2_queue_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src:." "${PYTHON}" - \
  "${ROOT}" "${STAGED}" "${LIVE}" "${STAGE}" <<'PY'
from hashlib import sha256
from pathlib import Path
import json
import os
import sys

from masld_bench.artifacts import freeze_tree
from scripts.gpu_bundle_dispatcher import load_item
from scripts.gpu_dependency_dispatcher import validate_item_route

root = Path(sys.argv[1]).resolve(strict=True)
staged = Path(sys.argv[2]).resolve(strict=True)
live = Path(sys.argv[3]).resolve(strict=False)
stage = Path(sys.argv[4]).resolve(strict=True)
item = validate_item_route(load_item(staged, root))
if (
    item["bundle_id"] != "model-training-333"
    or item["priority"] != 100
    or item["enabled"] is not True
    or item["family"] != "single_cell_chromatin"
    or item["logical_tasks"] != 10
    or item["_ready"] is not True
    or item["_job_name"] != "model-training-333"
    or item["_hours"] != 24.0
    or item["exports"]
    != {
        "PEAKVI_CROSS_COHORT_AUTHORITY": "peakvi_cross_cohort_v1",
        "PEAKVI_CROSS_COHORT_CONFIG": str(
            root
            / "config/campaigns/peakvi_gse296875_to_gse281367_fit_predict_revision2_20260825.json"
        ),
        "PEAKVI_CROSS_COHORT_CONFIG_SHA256": "a8c949c33469970d8f07fe98fd1464c2622d1ebee3f6928fa6213d62c2d3c795",
    }
):
    raise SystemExit("PeakVI dependency-queue route or identity differs")
bindings = {
    "campaign": (
        root
        / "config/campaigns/peakvi_gse296875_to_gse281367_fit_predict_revision2_20260825.json",
        "a8c949c33469970d8f07fe98fd1464c2622d1ebee3f6928fa6213d62c2d3c795",
    ),
    "admission": (
        root
        / "executions/peakvi-cross-cohort-revision2-admission-21101192/ARTIFACTS.json",
        "1fe974ab7c86b56186c899a3feccc26f71ce7fb3e7514cf40ee46ce25361c57a",
    ),
    "wrapper": (
        root / "slurm/run_peakvi_cross_cohort_revision2_bundle.sbatch",
        "dedff962d7a0847cd795d17ff888d71063eb4fc52e1104c9b083a8eb4110f553",
    ),
    "incident": (
        root / "config/artifacts/incidents/peakvi_target_cell_hash_collision_21101008.json",
        "97c0d926c020290ca9342026fbd4569350c083b62ffc7949de2290fd8332261f",
    ),
}
for label, (path, expected) in bindings.items():
    if sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit(f"PeakVI queue binding differs: {label}")
os.link(staged, live, follow_symlinks=False)
descriptor = os.open(live.parent, os.O_RDONLY)
try:
    os.fsync(descriptor)
finally:
    os.close(descriptor)
if sha256(live.read_bytes()).hexdigest() != item["_queue_sha256"]:
    raise SystemExit("published PeakVI queue hash differs")
receipt = {
    "schema_version": "masld-bench-peakvi-cross-cohort-queue-admission-v1",
    "status": "passed_and_published",
    "bundle_id": item["bundle_id"],
    "queue_sha256": item["_queue_sha256"],
    "wrapper_sha256": item["wrapper_sha256"],
    "priority": item["priority"],
    "logical_tasks": item["logical_tasks"],
    "requested_gpu_hours_ceiling": item["_hours"],
    "running_concurrency_policy": "four_afterany_dependency_lanes",
    "pending_job_limit": None,
    "account": "nslab",
    "qos": "nslab",
    "partition": "gpu",
    "gres": "gpu:l40s:1",
    "array_used": False,
    "innovation_used": False,
    "outcomes_read": False,
    "metrics_calculated": False,
    "scorer_invoked": False,
}
(stage / "receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(
    stage,
    {
        "artifact_class": "peakvi_cross_cohort_gpu_queue_admission",
        "bundle_id": item["bundle_id"],
        "queue_sha256": item["_queue_sha256"],
        "outcomes_read": False,
        "metrics_calculated": False,
        "status": "passed_and_published",
    },
)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
