#!/usr/bin/env bash
#SBATCH --job-name=model-check-289
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PLAN="${ROOT}/executions/sequence-five-seed-rectangle-admission-21098556"
TARGET="${ROOT}/executions/sequence-five-seed-queue-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite queue admission artifact.\n' >&2
  exit 2
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 78d785157d8b11b38d18f99362b98bdec2a894e6455f42ac54cf681958da0b28 \
  "${PLAN}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 1bad33ca17916aa1fde2683044b85fb67fb8b107fc536368ecc1cee1b20aa78c \
  "${ROOT}/scripts/materialize_sequence_task_native_five_seed_queue.py" | sha256sum --check --strict
printf '%s  %s\n' 7cc5201dca83652fc56dab912514ec9ca6d85f225d04ce65a8e4ff29dd6ed04b \
  "${ROOT}/scripts/gpu_bundle_dispatcher.py" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python \
  "${ROOT}/scripts/materialize_sequence_task_native_five_seed_queue.py" \
  --root "${ROOT}" \
  --plan-root "${PLAN}" \
  --output "${STAGE}/queue_items" \
  --priority-start 70
PYTHONPATH="${ROOT}/src:${ROOT}/scripts" python - "${ROOT}" "${STAGE}/queue_items" <<'PY'
from pathlib import Path
import json
import sys

from gpu_bundle_dispatcher import load_item

root = Path(sys.argv[1]).resolve(strict=True)
queue = Path(sys.argv[2]).resolve(strict=True)
paths = sorted(path for path in queue.glob("*.json") if path.name != "manifest.json")
items = [load_item(path, root) for path in paths]
manifest = json.loads((queue / "manifest.json").read_text(encoding="utf-8"))
if (
    len(items) != 30
    or manifest["queue_items"] != 30
    or manifest["priority_range"] != [70, 99]
    or any(not item["_ready"] for item in items)
    or any(item["family"] != "model_training" for item in items)
    or sum(item["logical_tasks"] for item in items) != 438
    or any(item["exports"]["MASLD_GPU_DISPATCHER_AUTHORITY"] != "sequence_five_seed_rectangle_v1" for item in items)
    or any("innovation" in json.dumps(item).lower() for item in items)
):
    raise SystemExit("materialized queue invariant failed")
print(json.dumps({
    "status": "pass",
    "queue_items": len(items),
    "logical_fits": sum(item["logical_tasks"] for item in items),
    "priority_range": manifest["priority_range"],
    "all_ready": all(item["_ready"] for item in items),
}, sort_keys=True))
PY
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${PLAN}/ARTIFACTS.json" \
  "${ROOT}/scripts/materialize_sequence_task_native_five_seed_queue.py" \
  "${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "${ROOT}/slurm/admit_sequence_task_native_five_seed_queue_cpu.sbatch" \
  > "${STAGE}/validation/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = json.loads((root / "queue_items/manifest.json").read_text(encoding="utf-8"))
freeze_tree(root, {
    "artifact_class": "sequence_task_native_five_seed_queue_admission",
    "status": "passed",
    "queue_items": manifest["queue_items"],
    "priority_range": manifest["priority_range"],
    "plan_artifacts_sha256": manifest["plan_artifacts_sha256"],
    "submission_authority": "dispatcher_only_no_manual_gpu_sbatch",
    "qos": "nslab",
    "gpu_squeue_ceiling": 5,
    "manual_gpu_sbatch_allowed": False,
    "prediction_or_metric_values_read": False,
    "outcome_paths_available": False,
    "sealed_paths_available": False,
})
PY
module purge
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
