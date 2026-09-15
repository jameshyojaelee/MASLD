#!/usr/bin/env bash
#SBATCH --job-name=validate-model-work-203
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
QUEUE="${ROOT}/config/campaigns/gpu_bundle_queue/027-model-work-203-retry1.json"
TARGET="${ROOT}/executions/model-work-203-queue-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}/src" python - "${QUEUE}" "${STAGE}/validation/receipt.json" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

queue_path, receipt_path = map(Path, sys.argv[1:])
root = Path.cwd().resolve(strict=True)
queue = json.loads(queue_path.read_text())
expected = {
    "schema_version": "masld-bench-gpu-bundle-queue-item-v1",
    "bundle_id": "model-work-203-retry1",
    "priority": 27,
    "enabled": True,
    "wrapper_path": "slurm/run_model_work_203.sbatch",
    "wrapper_sha256": "97fdcbeb58399535e379dd53b648d8b48551704be2c562ab6e2567d8b4aa705e",
    "exports": {},
    "logical_tasks": 1033,
    "family": "sequence_regulatory",
}
for key, value in expected.items():
    if queue.get(key) != value:
        raise SystemExit(f"queue contract differs: {key}")
wrapper = root / queue["wrapper_path"]
if sha256(wrapper.read_bytes()).hexdigest() != queue["wrapper_sha256"]:
    raise SystemExit("GPU wrapper hash differs")
text = wrapper.read_text()
for required in (
    "#SBATCH --job-name=model-work-203",
    "#SBATCH --partition=gpu",
    "#SBATCH --gres=gpu:l40s:1",
    "#SBATCH --qos=nslab",
    "--containall --net --network none",
    '"${ROOT}:${ROOT}:ro"',
    '"${STAGE}:${STAGE}:rw"',
    '"${CACHE}:${CACHE}:rw"',
    '--pwd "${ROOT}"',
    '--env PYTHONPATH="${ROOT}"',
    '"${FIXTURE_AUDIT}"',
    '"${PREFLIGHT}"',
):
    if text.count(required) != 1:
        raise SystemExit(f"GPU routing differs: {required}")
for forbidden in ("--qos=innovation", "PLACEHOLDER_", "masld-enformer"):
    if forbidden in text:
        raise SystemExit(f"forbidden wrapper value: {forbidden}")
for relative in queue.get("required_paths", []):
    artifact = (root / relative).resolve(strict=True)
    if artifact.name != "ARTIFACTS.json" or root not in artifact.parents:
        raise SystemExit("required artifact path differs")
    verify_frozen_tree(artifact.parent)
receipt = {
    "schema_version": "masld-bench-model-work-203-queue-validation-v1",
    "status": "pass_enabled_queue_item",
    "bundle_id": queue["bundle_id"],
    "enabled": queue["enabled"],
    "logical_tasks": queue["logical_tasks"],
    "wrapper_sha256": queue["wrapper_sha256"],
    "queue_sha256": sha256(queue_path.read_bytes()).hexdigest(),
    "gpu_submitted": False,
    "outcomes_read": False,
}
receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
PY
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum config/campaigns/gpu_bundle_queue/027-model-work-203-retry1.json \
  slurm/run_model_work_203.sbatch slurm/validate_model_work_203_queue.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/receipt.json").read_text())
if receipt["status"] != "pass_enabled_queue_item" or not receipt["enabled"] or receipt["gpu_submitted"]:
    raise SystemExit("queue validation receipt differs")
freeze_tree(root, {
    "artifact_class": "enformer_gse281364_gpu_queue_validation",
    "bundle_id": "model-work-203-retry1",
    "enabled": True,
    "gpu_submitted": False,
    "outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
