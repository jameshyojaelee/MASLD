#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-216
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-216-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-216-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BUNDLE=${ROOT}/executions/dispatch-bundles/cobolt-5seed-25run-v9-generic6
BUNDLE_SHA=3015e5179affe11f5b2edb16b4197d75cdcf7adb27f0269388985b24c06597ad
WRAPPER=${BUNDLE}/model_work_212.sbatch
WRAPPER_SHA=ba30e7d376e2a42ee65b6b0a8958f675e03090faa87dac9fbccd10abe0aa273a
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58
CANDIDATE_SHA=ac24e1a8da4fde1c35aafa1370659155eded36811099b40f3ee39a4fc2866862
STABILITY=${ROOT}/executions/cobolt-v9-frozen-candidate-stability-21087850
STABILITY_SHA=b8b6003576becb93b530c20465d9022ee334bae5bbdcb660ba340836d0831233
ITEM=${ROOT}/config/campaigns/gpu_bundle_queue/057-cobolt-5seed-25run-v9-generic6.json
TARGET=${ROOT}/executions/cobolt-v9-generic6-gpu-queue-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite Cobolt generic6 queue validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
for binding in \
  "${BUNDLE_SHA}  ${BUNDLE}/ARTIFACTS.json" \
  "${WRAPPER_SHA}  ${WRAPPER}" \
  "${CANDIDATE_SHA}  ${CANDIDATE}/ARTIFACTS.json" \
  "${STABILITY_SHA}  ${STABILITY}/ARTIFACTS.json"; do
  printf '%s\n' "${binding}" | sha256sum --check --strict
done
export PYTHONPATH=${ROOT}:${ROOT}/src
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
for tree in (bundle, candidate, stability):
    verify_frozen_tree(tree)
manifest = json.loads((bundle / "bundle_manifest.json").read_text())
plan = load_frozen_plan(candidate)
stability_receipt = json.loads(
    (stability / "resource_stability_receipt.json").read_text()
)
wrapper = (bundle / "model_work_212.sbatch").read_text()
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
    manifest.get("bundle_id") != "cobolt-5seed-25run-v9-generic6"
    or not item["_ready"]
    or item["bundle_id"] != manifest.get("bundle_id")
    or item["wrapper_sha256"] != sha256_file(bundle / "model_work_212.sbatch")
    or item["exports"] != {}
    or item["logical_tasks"] != 25
    or manifest.get("resources") != expected_resources
    or manifest.get("candidate_manifest_sha256") != sha256_file(candidate / "ARTIFACTS.json")
    or manifest.get("plan_sha256") != plan.get("plan_sha256")
    or len(manifest.get("runs", [])) != 25
    or stability_receipt.get("candidate_artifacts_sha256")
    != sha256_file(candidate / "ARTIFACTS.json")
    or stability_receipt.get("status") != "pass"
    or stability_receipt.get("outcomes_opened") is not False
    or "#SBATCH --job-name=model-work-212" not in wrapper
    or "#SBATCH --partition=gpu" not in wrapper
    or "#SBATCH --account=nslab" not in wrapper
    or "#SBATCH --qos=nslab" not in wrapper
    or "innovation" in wrapper.lower()
    or "masld" in next(
        line.lower() for line in wrapper.splitlines() if "--job-name=" in line
    )
):
    raise SystemExit("Cobolt generic6 GPU queue contract differs")
write_json_exclusive(stage / "queue_validation_receipt.json", {
    "schema_version": "masld-bench-cobolt-v9-generic6-gpu-queue-validation-v1",
    "status": "pass",
    "bundle_id": manifest["bundle_id"],
    "queue_item_sha256": item["_queue_sha256"],
    "bundle_artifacts_sha256": sha256_file(bundle / "ARTIFACTS.json"),
    "wrapper_sha256": sha256_file(bundle / "model_work_212.sbatch"),
    "candidate_artifacts_sha256": sha256_file(candidate / "ARTIFACTS.json"),
    "stability_artifacts_sha256": sha256_file(stability / "ARTIFACTS.json"),
    "logical_tasks": 25,
    "gpu_jobs_submitted": 0,
    "outcomes_opened": False,
})
freeze_tree(stage, {
    "artifact_class": "cobolt_v9_generic6_gpu_queue_validation",
    "bundle_id": manifest["bundle_id"],
    "gpu_jobs_submitted": 0,
    "outcomes_opened": False,
    "status": "pass",
})
verify_frozen_tree(stage)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
