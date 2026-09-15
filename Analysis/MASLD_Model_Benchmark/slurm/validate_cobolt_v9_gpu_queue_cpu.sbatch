#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-213
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-213-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-213-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ITEM=${ROOT}/config/campaigns/gpu_bundle_queue_staging/056-cobolt-5seed-25run-v9.json
BUNDLE=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v9-generic5
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--7031e0dfbf76bc4d
STABILITY=${ROOT}/executions/cobolt-v9-resource-stability-21083320
TARGET=${ROOT}/executions/cobolt-v9-gpu-queue-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Cobolt v9 queue validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python -m py_compile scripts/gpu_bundle_dispatcher.py scripts/build_cobolt_v9_gpu_bundle.py
python - "${ROOT}" "${ITEM}" "${BUNDLE}" "${CANDIDATE}" "${STABILITY}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file
from masld_bench.planner import load_frozen_plan
from scripts.gpu_bundle_dispatcher import load_item

root, item_path, bundle, candidate, stability, stage = map(
    lambda value: Path(value).resolve(strict=True), sys.argv[1:]
)
item = load_item(item_path, root)
if not item["_ready"] or item["bundle_id"] != "cobolt-5seed-25run-v9-generic5":
    raise SystemExit("Cobolt v9 staged queue item is not ready")
for tree in (bundle, candidate, stability):
    verify_frozen_tree(tree)
manifest = json.loads((bundle / "bundle_manifest.json").read_text(encoding="utf-8"))
plan = load_frozen_plan(candidate)
stability_receipt = json.loads(
    (stability / "resource_stability_receipt.json").read_text(encoding="utf-8")
)
expected_resources = {
    "account": "nslab",
    "partition": "gpu",
    "qos": "nslab",
    "gpu_type": "l40s",
    "gpus": 1,
    "cpus": 4,
    "memory_gb": 64,
    "walltime": "06:00:00",
}
if (
    manifest.get("resources") != expected_resources
    or manifest.get("candidate_manifest_sha256") != sha256_file(candidate / "ARTIFACTS.json")
    or manifest.get("plan_sha256") != plan.get("plan_sha256")
    or len(manifest.get("runs", [])) != 25
    or stability_receipt.get("status") != "pass"
    or stability_receipt.get("outcomes_opened") is not False
):
    raise SystemExit("Cobolt v9 queue bindings differ")
write_json_exclusive(stage / "queue_validation_receipt.json", {
    "schema_version": "masld-bench-cobolt-v9-gpu-queue-validation-v1",
    "status": "pass",
    "bundle_id": item["bundle_id"],
    "queue_item_sha256": item["_queue_sha256"],
    "bundle_artifacts_sha256": sha256_file(bundle / "ARTIFACTS.json"),
    "candidate_artifacts_sha256": sha256_file(candidate / "ARTIFACTS.json"),
    "stability_artifacts_sha256": sha256_file(stability / "ARTIFACTS.json"),
    "logical_tasks": 25,
    "gpu_jobs_submitted": 0,
    "outcomes_opened": False,
})
freeze_tree(stage, {
    "artifact_class": "cobolt_v9_gpu_queue_validation",
    "bundle_id": item["bundle_id"],
    "status": "pass",
    "gpu_jobs_submitted": 0,
    "outcomes_opened": False,
})
verify_frozen_tree(stage)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
