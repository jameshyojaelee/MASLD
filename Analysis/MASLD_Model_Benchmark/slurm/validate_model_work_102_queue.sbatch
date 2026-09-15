#!/usr/bin/env bash
#SBATCH --job-name=model-work-102-queue-check
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
QUEUE="${ROOT}/config/campaigns/gpu_bundle_queue/024-model-work-102.json"
TARGET="${ROOT}/executions/model-work-102-queue-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite queue validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' \
  ba8e3ee9183bdb7723ee7a454d3b2052814cc6b781833848fd8993d91e62e784 \
  "${QUEUE}" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ROOT}" "${QUEUE}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root, queue, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(root))
from scripts.gpu_bundle_dispatcher import load_item

item = load_item(queue, root.resolve(strict=True))
if (
    item["bundle_id"] != "model-work-102"
    or item["enabled"]
    or not item["_ready"]
    or item["wrapper_sha256"] != "363de52e6ceda82d9908d3a24fed47a18a164613eb16126c24ffbef0e564a225"
    or item["logical_tasks"] != 5
    or item["family"] != "local_variant_effect"
):
    raise SystemExit("disabled model-work-102 queue contract differs")
receipt = {
    "schema_version": "masld-bench-model-work-102-queue-validation-v1",
    "status": "pass_disabled_ready_queue_item",
    "bundle_id": item["bundle_id"],
    "enabled": item["enabled"],
    "ready": item["_ready"],
    "wrapper_sha256": item["wrapper_sha256"],
    "queue_sha256": item["_queue_sha256"],
    "logical_tasks": item["logical_tasks"],
    "gpu_submission_performed": False,
    "outcomes_read": False,
}
(output / "validation/queue_receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(receipt, sort_keys=True))
PY
sha256sum \
  "${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "${QUEUE}" \
  "${ROOT}/slurm/run_model_work_102.sbatch" \
  "${ROOT}/slurm/validate_model_work_102_queue.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/queue_receipt.json").read_text())
if (
    receipt["status"] != "pass_disabled_ready_queue_item"
    or receipt["enabled"]
    or not receipt["ready"]
    or receipt["gpu_submission_performed"]
    or receipt["outcomes_read"]
):
    raise SystemExit("queue validation receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "disabled_gpu_queue_validation",
        "bundle_id": "model-work-102",
        "enabled": False,
        "ready": True,
        "gpu_submission_performed": False,
        "outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-work-102 disabled queue validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
