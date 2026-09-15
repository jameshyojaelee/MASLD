#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-queue-205
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-queue-205-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-queue-205-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
QUEUE=${ROOT}/config/campaigns/gpu_bundle_queue/052-cobolt-5seed-25run-v8.json
TARGET=${ROOT}/executions/cobolt-v8-gpu-queue-validation-v1
STAGE=${TARGET}.staging-${SLURM_JOB_ID}
FAILED=${TARGET}.failed-${SLURM_JOB_ID}
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Cobolt v8 queue validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}:${ROOT}/src" python - \
  "${QUEUE}" "${STAGE}/validation/queue_receipt.json" <<'PY'
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

from scripts.gpu_bundle_dispatcher import load_item
from masld_bench.artifacts import verify_frozen_tree

queue_path = Path(sys.argv[1]).resolve(strict=True)
receipt_path = Path(sys.argv[2])
root = Path.cwd().resolve(strict=True)
item = load_item(queue_path, root)

expected_exports = {
    "COBOLT_V8_BUNDLE_SHA256": "b5a41bbd21458d5b60dbda415203c2cf07febee8e747801160cbef7b63941a7e",
    "COBOLT_V8_CANDIDATE_SHA256": "d6fd234073768e7f2b96389aa420037df5e375c16e3ddeb742b9b101a030f8c0",
    "COBOLT_V8_CPU_VERIFIER_SHA256": "da36ab47966f93f52deb9372838e47cbe5ca89bc5463b6a4bc790aef4410279d",
    "COBOLT_V8_EVALUATOR_SOURCE_SHA256": "cce3a0be9b50b036d6edc971fefe09af7e90aa12cc957ae18b45bf9c63251e7e",
}
expected_required = {
    "executions/dispatch-bundles/cobolt-5seed-25run-v8-generic5/ARTIFACTS.json": expected_exports["COBOLT_V8_BUNDLE_SHA256"],
    "candidates/v1-rna-atac-cobolt-5seed-smoke-v8--720da062fc1b50ab/ARTIFACTS.json": expected_exports["COBOLT_V8_CANDIDATE_SHA256"],
    "executions/cobolt-v8-cpu-baseline-verification-21078928/ARTIFACTS.json": expected_exports["COBOLT_V8_CPU_VERIFIER_SHA256"],
    "executions/rna-atac-cobolt-5seed-v8-evaluator-source-21077948/ARTIFACTS.json": expected_exports["COBOLT_V8_EVALUATOR_SOURCE_SHA256"],
    "executions/cobolt-v8-gpu-queue-validation-v1/ARTIFACTS.json": None,
}
if (
    item["bundle_id"] != "cobolt-5seed-25run-v8-generic5"
    or item["priority"] != 52
    or not item["enabled"]
    or item["_ready"]
    or item["logical_tasks"] != 25
    or item["family"] != "observed_multiome"
    or item["wrapper_sha256"] != "8440b47e1de682a26f83d93adc5cbc6825a02638723906d0463a84faef062beb"
    or item["exports"] != expected_exports
    or item["required_paths"] != list(expected_required)
):
    raise SystemExit("Cobolt v8 pre-activation queue contract differs")

for relative, expected in expected_required.items():
    if expected is None:
        continue
    artifact = (root / relative).resolve(strict=True)
    if sha256(artifact.read_bytes()).hexdigest() != expected:
        raise SystemExit(f"required artifact hash differs: {relative}")
    verify_frozen_tree(artifact.parent)

wrapper = item["_wrapper"].read_text(encoding="utf-8")
for required in (
    "#SBATCH --job-name=model-work-205",
    "#SBATCH --account=nslab",
    "#SBATCH --partition=gpu",
    "#SBATCH --qos=nslab",
    "#SBATCH --cpus-per-task=4",
    "#SBATCH --mem=64G",
    "#SBATCH --time=06:00:00",
    "#SBATCH --gres=gpu:l40s:1",
):
    if wrapper.count(required) != 1:
        raise SystemExit(f"Cobolt v8 wrapper differs: {required}")
for forbidden in ("#SBATCH --array", "--qos=innovation", "#SBATCH --job-name=masld"):
    if forbidden in wrapper:
        raise SystemExit(f"forbidden Cobolt v8 wrapper value: {forbidden}")

bundle_root = root / "executions/dispatch-bundles/cobolt-5seed-25run-v8-generic5"
manifest = json.loads((bundle_root / "bundle_manifest.json").read_text(encoding="utf-8"))
runs = manifest["runs"]
if (
    manifest["candidate_manifest_sha256"] != expected_exports["COBOLT_V8_CANDIDATE_SHA256"]
    or manifest["resources"] != {
        "account": "nslab", "partition": "gpu", "qos": "nslab",
        "gpu_type": "l40s", "gpus": 1, "cpus": 4,
        "memory_gb": 64, "walltime": "06:00:00",
    }
    or len(runs) != 25
    or len({run["run_id"] for run in runs}) != 25
    or {(run["seed"], run["fold"]) for run in runs}
        != {(seed, fold) for seed in (1103, 2107, 3109, 4111, 5113) for fold in range(5)}
):
    raise SystemExit("Cobolt v8 bundle census or resources differ")
continuation = manifest["continuation_policy"]
if (
    continuation.get("mid_run_resume") is not False
    or continuation.get("interrupted_run_recovery") != "new_campaign_revision_required"
    or continuation.get("bundle_retry_identity")
        != "distinct_immutable_bundle_and_queue_item_required"
):
    raise SystemExit("Cobolt v8 retry limitation differs")

states = subprocess.run(
    ["squeue", "-u", os.environ["USER"], "-p", "gpu", "-h", "-o", "%T"],
    check=True,
    capture_output=True,
    text=True,
).stdout.splitlines()
receipt = {
    "schema_version": "masld-bench-cobolt-v8-gpu-queue-validation-v1",
    "status": "pass_activation_artifact",
    "bundle_id": item["bundle_id"],
    "enabled": True,
    "ready_before_validation_artifact": False,
    "ready_after_validation_artifact": True,
    "logical_tasks": 25,
    "queue_sha256": item["_queue_sha256"],
    "wrapper_sha256": item["wrapper_sha256"],
    "required_artifact_sha256": expected_exports,
    "gpu_states_before_activation": {
        "running": sum(state == "RUNNING" for state in states),
        "pending": sum(state == "PENDING" for state in states),
        "total": len(states),
    },
    "gpu_submission_performed_by_validator": False,
    "mid_run_resume": False,
    "interrupted_run_recovery": "new_campaign_revision_required",
    "outcomes_read": False,
}
receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${QUEUE}" \
  "${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "${ROOT}/slurm/validate_cobolt_v8_gpu_queue.sbatch" \
  "${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v8-generic5/model_work_205.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/queue_receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_activation_artifact"
    or not receipt["enabled"]
    or receipt["ready_before_validation_artifact"]
    or not receipt["ready_after_validation_artifact"]
    or receipt["gpu_submission_performed_by_validator"]
    or receipt["outcomes_read"]
):
    raise SystemExit("Cobolt v8 queue validation receipt differs")
freeze_tree(root, {
    "artifact_class": "cobolt_v8_gpu_queue_activation_validation",
    "bundle_id": receipt["bundle_id"],
    "gpu_submission_performed_by_validator": False,
    "mid_run_resume": False,
    "outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Cobolt v8 GPU queue validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
