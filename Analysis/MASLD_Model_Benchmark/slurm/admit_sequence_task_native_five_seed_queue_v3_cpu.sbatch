#!/usr/bin/env bash
#SBATCH --job-name=model-check-292
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
CORRECTION="${ROOT}/executions/sequence-five-seed-execution-correction-21098951"
TARGET="${ROOT}/executions/sequence-five-seed-queue-admission-v3-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite corrected queue admission artifact.\n' >&2
  exit 2
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 78d785157d8b11b38d18f99362b98bdec2a894e6455f42ac54cf681958da0b28 \
  "${PLAN}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 5f3c122972f86a16e615c52c98e1f20af15fc6dcbf2c83b4b0f53730b0a68176 \
  "${CORRECTION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' d0d2dd5264a31ca6cfa7908aff70957ddc4e432f30ef773dd1d3318bd968e52b \
  "${ROOT}/scripts/materialize_sequence_task_native_five_seed_queue_v2.py" | sha256sum --check --strict
printf '%s  %s\n' 7cc5201dca83652fc56dab912514ec9ca6d85f225d04ce65a8e4ff29dd6ed04b \
  "${ROOT}/scripts/gpu_bundle_dispatcher.py" | sha256sum --check --strict
printf '%s  %s\n' 7cfd70d94ad9fb57c4ae6337ed3696e8c4a381bc5d5a96b6800e8ef744c13e3d \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python \
  "${ROOT}/scripts/materialize_sequence_task_native_five_seed_queue_v2.py" \
  --root "${ROOT}" \
  --plan-root "${PLAN}" \
  --correction-root "${CORRECTION}" \
  --output "${STAGE}/queue_items" \
  --priority-start 70
PYTHONPATH="${ROOT}/src:${ROOT}/scripts" python - \
  "${ROOT}" "${STAGE}/queue_items" "${PLAN}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from gpu_bundle_dispatcher import load_item

root = Path(sys.argv[1]).resolve(strict=True)
generated = Path(sys.argv[2]).resolve(strict=True)
plan = Path(sys.argv[3]).resolve(strict=True)
live = root / "config/campaigns/gpu_bundle_queue"
claims = root / "executions/gpu-bundle-dispatch-state/claims"
paths = sorted(path for path in generated.glob("*.json") if path.name != "manifest.json")
items = [load_item(path, root) for path in paths]
manifest = json.loads((generated / "manifest.json").read_text(encoding="utf-8"))
live_items = [json.loads(path.read_text(encoding="utf-8")) for path in live.glob("*.json")]
live_ids = {item["bundle_id"] for item in live_items}
v1 = [item for item in live_items if item["bundle_id"] in {f"model-training-{n}" for n in range(701, 731)}]
generated_ids = {item["bundle_id"] for item in items}
task_ids = {item["exports"]["RECTANGLE_BUNDLE_ID"] for item in items}
with (plan / "plan/bundle_allocations.tsv").open(encoding="utf-8", newline="") as handle:
    expected_task_ids = {row["bundle_id"] for row in csv.DictReader(handle, delimiter="\t")}
if (
    len(items) != 30
    or len(generated_ids) != 30
    or generated_ids & live_ids
    or any((claims / f"{bundle_id}.json").exists() for bundle_id in generated_ids)
    or len(v1) != 30
    or any(item["enabled"] for item in v1)
    or any((claims / f"{item['bundle_id']}.json").exists() for item in v1)
    or task_ids != expected_task_ids
    or manifest["logical_fits"] != 438
    or manifest["priority_range"] != [70, 99]
    or any(not item["_ready"] for item in items)
    or any(item["family"] != "model_training" for item in items)
    or sum(item["logical_tasks"] for item in items) != 438
    or any("innovation" in path.read_text(encoding="utf-8").lower() for path in paths)
):
    raise SystemExit("corrected dispatcher identity, collision, or no-v1-claim invariant failed")
print(json.dumps({
    "status": "pass",
    "queue_items": len(items),
    "logical_fits": sum(item["logical_tasks"] for item in items),
    "dispatcher_identity_noncollision": True,
    "v1_disabled_unclaimed": True,
    "all_ready": all(item["_ready"] for item in items),
}, sort_keys=True))
PY
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${PLAN}/ARTIFACTS.json" \
  "${CORRECTION}/ARTIFACTS.json" \
  "${ROOT}/scripts/materialize_sequence_task_native_five_seed_queue_v2.py" \
  "${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch" \
  "${ROOT}/slurm/admit_sequence_task_native_five_seed_queue_v3_cpu.sbatch" \
  > "${STAGE}/validation/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = json.loads((root / "queue_items/manifest.json").read_text(encoding="utf-8"))
freeze_tree(root, {
    "artifact_class": "sequence_task_native_five_seed_corrected_queue_admission",
    "status": "passed",
    "queue_items": manifest["queue_items"],
    "logical_fits": manifest["logical_fits"],
    "priority_range": manifest["priority_range"],
    "dispatcher_identity_range": manifest["dispatcher_identity_range"],
    "task_identity_range": manifest["task_identity_range"],
    "plan_artifacts_sha256": manifest["plan_artifacts_sha256"],
    "correction_artifacts_sha256": manifest["correction_artifacts_sha256"],
    "v1_queue_items_disabled_unclaimed": True,
    "dispatcher_identity_noncollision": True,
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
