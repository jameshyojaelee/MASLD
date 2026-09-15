#!/usr/bin/env bash
#SBATCH --job-name=model-check-291
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
TARGET="${ROOT}/executions/sequence-five-seed-execution-correction-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite correction admission artifact.\n' >&2
  exit 2
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 78d785157d8b11b38d18f99362b98bdec2a894e6455f42ac54cf681958da0b28 \
  "${PLAN}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 1f9bcc4e1a91628753b2f4ba144ca5f69bf317ce0be2c694ede51773d823bba8 \
  "${ROOT}/config/campaigns/sequence_task_native_five_seed_execution_correction_20260825.json" | sha256sum --check --strict
printf '%s  %s\n' 4da56214214fd2e876c9b4dd9d9d49746782f3a972fe62bf2bc7bd12208a0582 \
  "${ROOT}/scripts/run_sequence_task_native_five_seed_bundle.py" | sha256sum --check --strict
printf '%s  %s\n' a57f8b222da3e2c6179d9300eb58a1186af91a4b87274281a239efbe7a4af6c7 \
  "${ROOT}/scripts/run_sequence_task_native_five_seed_bundle_v2.py" | sha256sum --check --strict
printf '%s  %s\n' d55fd67898661a8f652ee6ee74d5fcc07b01b4b03d96d08d1a79c9276f9295b2 \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle.sbatch" | sha256sum --check --strict
printf '%s  %s\n' 7cfd70d94ad9fb57c4ae6337ed3696e8c4a381bc5d5a96b6800e8ef744c13e3d \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch" | sha256sum --check --strict
printf '%s  %s\n' 7061307461994518585f25470e2241eb44b5fd116c6b994db788b1018b85ed69 \
  "${ROOT}/tests/unit/test_sequence_task_native_five_seed_execution_correction.py" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python -m unittest -v \
  tests.unit.test_sequence_task_native_five_seed_execution_correction \
  > "${STAGE}/validation/unit.stdout.txt" \
  2> "${STAGE}/validation/unit.stderr.txt"
python - "${ROOT}" "${STAGE}/validation/queue-disposition.json" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
queue = root / "config/campaigns/gpu_bundle_queue"
claims = root / "executions/gpu-bundle-dispatch-state/claims"
paths = sorted(queue.glob("0[7-9][0-9]-model-training-7*.json"))
items = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
claimed = [item["bundle_id"] for item in items if (claims / f"{item['bundle_id']}.json").exists()]
if len(items) != 30 or claimed or any(item["enabled"] for item in items):
    raise SystemExit("v1 queue was not disabled before correction admission")
Path(sys.argv[2]).write_text(json.dumps({
    "status": "pass",
    "v1_queue_items": len(items),
    "v1_claimed": claimed,
    "v1_submitted": [],
    "v1_enabled": 0,
}, sort_keys=True, indent=2) + "\n", encoding="utf-8")
PY
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${PLAN}/ARTIFACTS.json" \
  "${ROOT}/config/campaigns/sequence_task_native_five_seed_execution_correction_20260825.json" \
  "${ROOT}/scripts/run_sequence_task_native_five_seed_bundle.py" \
  "${ROOT}/scripts/run_sequence_task_native_five_seed_bundle_v2.py" \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle.sbatch" \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch" \
  "${ROOT}/tests/unit/test_sequence_task_native_five_seed_execution_correction.py" \
  > "${STAGE}/validation/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
disposition = json.loads((root / "validation/queue-disposition.json").read_text(encoding="utf-8"))
freeze_tree(root, {
    "artifact_class": "sequence_task_native_five_seed_execution_correction",
    "status": "passed",
    "rectangle_plan_artifacts_sha256": "78d785157d8b11b38d18f99362b98bdec2a894e6455f42ac54cf681958da0b28",
    "v1_queue_items_disabled": disposition["v1_queue_items"],
    "v1_queue_items_claimed": len(disposition["v1_claimed"]),
    "v1_queue_items_submitted": len(disposition["v1_submitted"]),
    "controller_diff": "Counter_import_plus_testable_terminal_summary_extraction",
    "scientific_task_contract_changed": False,
    "prediction_or_metric_values_read": False,
    "outcome_paths_available": False,
    "sealed_paths_available": False,
    "manual_gpu_submission_allowed": False,
    "qos": "nslab",
})
PY
module purge
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
