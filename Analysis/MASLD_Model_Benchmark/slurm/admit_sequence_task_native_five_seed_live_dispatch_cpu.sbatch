#!/usr/bin/env bash
#SBATCH --job-name=model-check-294
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
QUEUE_ADMISSION="${ROOT}/executions/sequence-five-seed-queue-admission-v3-21098963"
CONTRACT="${ROOT}/config/campaigns/sequence_task_native_five_seed_live_dispatch_reconciliation_20260825.json"
TARGET="${ROOT}/executions/sequence-five-seed-live-dispatch-readiness-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite five-seed live-dispatch readiness artifact.\n' >&2
  exit 2
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' 78d785157d8b11b38d18f99362b98bdec2a894e6455f42ac54cf681958da0b28 \
  "${PLAN}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 5f3c122972f86a16e615c52c98e1f20af15fc6dcbf2c83b4b0f53730b0a68176 \
  "${CORRECTION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ceca0985719a13b7588dbc5156401f30b611906338aecc0bfaa07ce0eba5a646 \
  "${QUEUE_ADMISSION}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH="${ROOT}/src:${ROOT}/scripts" python -m unittest -v \
  tests.unit.test_sequence_task_native_five_seed_execution_correction \
  tests.unit.test_sequence_task_native_five_seed_live_dispatch \
  > "${STAGE}/validation/unit.stdout.txt" \
  2> "${STAGE}/validation/unit.stderr.txt"
PYTHONPATH="${ROOT}/src" python \
  "${ROOT}/scripts/audit_sequence_task_native_five_seed_live_dispatch.py" \
  --root "${ROOT}" \
  --contract "${CONTRACT}" \
  --plan-root "${PLAN}" \
  --correction-root "${CORRECTION}" \
  --queue-admission-root "${QUEUE_ADMISSION}" \
  --output "${STAGE}/reconciliation" \
  --user "${USER}" \
  > "${STAGE}/validation/audit.stdout.json"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONTRACT}" \
  "${ROOT}/scripts/audit_sequence_task_native_five_seed_live_dispatch.py" \
  "${ROOT}/tests/unit/test_sequence_task_native_five_seed_live_dispatch.py" \
  "${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "${ROOT}/scripts/run_sequence_task_native_five_seed_bundle_v2.py" \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch" \
  "${ROOT}/slurm/admit_sequence_task_native_five_seed_live_dispatch_cpu.sbatch" \
  > "${STAGE}/validation/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "reconciliation/summary.json").read_text(encoding="utf-8"))
if (
    summary["status"] != "pass"
    or summary["v1_enabled"] != 0
    or summary["v1_claimed"] != 0
    or summary["v1_ledger_submissions"] != 0
    or summary["v2_enabled_ready"] != 30
    or summary["logical_fits"] != 438
    or summary["central_dispatcher_running"] is not True
    or summary["manual_gpu_submission_performed"] is not False
    or summary["qos"] != "nslab"
    or summary["innovation_used"] is not False
    or summary["prediction_or_metric_values_read"] is not False
    or summary["outcomes_read"] is not False
    or summary["sealed_data_read"] is not False
):
    raise SystemExit("live-dispatch readiness receipt differs")
freeze_tree(root, {
    "artifact_class": "sequence_task_native_five_seed_live_dispatch_readiness",
    "status": "passed",
    "readiness": summary["readiness"],
    "v1_queue_items": summary["v1_queue_items"],
    "v1_enabled": summary["v1_enabled"],
    "v1_claimed": summary["v1_claimed"],
    "v1_ledger_submissions": summary["v1_ledger_submissions"],
    "v2_queue_items": summary["v2_queue_items"],
    "v2_enabled_ready": summary["v2_enabled_ready"],
    "v2_claimed": summary["v2_claimed"],
    "logical_fits": summary["logical_fits"],
    "controller_import_complete": summary["counter_import_complete"],
    "terminal_summary_paths_executed": True,
    "dispatcher_bundle_identities_distinct": True,
    "dispatcher_bundle_ids_globally_unique": True,
    "live_v2_byte_identical_to_frozen_admission": True,
    "central_dispatcher_running": True,
    "gpu_jobs": summary["gpu_jobs"],
    "gpu_running_or_other_active": summary["gpu_running_or_other_active"],
    "gpu_pending": summary["gpu_pending"],
    "gpu_job_ceiling": 5,
    "maximum_running": 4,
    "maximum_pending": 1,
    "submission_authority": "central_dispatcher_only",
    "manual_gpu_submission_performed": False,
    "account": "nslab",
    "qos": "nslab",
    "innovation_used": False,
    "prediction_or_metric_values_read": False,
    "outcomes_read": False,
    "sealed_data_read": False,
})
PY
module purge
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
