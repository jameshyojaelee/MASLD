#!/usr/bin/env bash
#SBATCH --job-name=model-check-279
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python3.11
QUEUE="${ROOT}/config/campaigns/gpu_bundle_queue/061-model-training-603.json"
TARGET="${ROOT}/executions/model-check-279-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-training-603 queue validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "3424ee597c8873bde03e00d6a7f5bc9a415933990c6381e8f8ebd769d327465a ${QUEUE}" \
  "ac152c23e45890400576421b3b4b47fa015c890b459eae5b6440da006807c606 ${ROOT}/slurm/run_model_training_603.sbatch" \
  "7c81f285a121b9d70297ad7406bb40e445dda6585478f5896de9c41827d2e894 ${ROOT}/config/borzoi_grelu_converted_port_full_window_probe.json" \
  "d9de47bceea5e32d5da1d70b0ac66155984eeca56c2c655c6f30056200969bb6 ${ROOT}/scripts/probe_borzoi_grelu_converted_port_full_window.py" \
  "d6edbe8c53f2caa7fce29cb4768b0f98dbdfc9f69a01bc5b9e78f9661b327627 ${ROOT}/tests/unit/test_probe_borzoi_grelu_converted_port_full_window.py" \
  "7cc5201dca83652fc56dab912514ec9ca6d85f225d04ce65a8e4ff29dd6ed04b ${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "06eb07080d4cae46eeb373fd309dadc971e149a74ce857adea1c3945543e901e ${ROOT}/tests/unit/test_gpu_bundle_dispatcher.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_probe_borzoi_grelu_converted_port_full_window \
  > "${STAGE}/validation/model-tests.stdout.txt" \
  2> "${STAGE}/validation/model-tests.stderr.txt"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}" python - "${QUEUE}" "${ROOT}" \
  > "${STAGE}/validation/queue.json" <<'PY'
import json
from pathlib import Path
import sys
from scripts.gpu_bundle_dispatcher import load_item

item = load_item(Path(sys.argv[1]), Path(sys.argv[2]))
if (
    item["bundle_id"] != "model-training-603"
    or item["priority"] != 61
    or item["enabled"]
    or not item["_ready"]
    or item["logical_tasks"] != 1
    or item["family"] != "long_range_regulatory"
):
    raise SystemExit("disabled model-training-603 queue item differs")
print(json.dumps({
    "bundle_id": item["bundle_id"], "priority": item["priority"],
    "enabled": item["enabled"], "ready": item["_ready"],
    "queue_sha256": item["_queue_sha256"],
    "wrapper_sha256": item["wrapper_sha256"],
}, indent=2, sort_keys=True))
PY
MASLD_GPU_DISPATCHER_AUTHORITY=sequence_gpu_dispatcher_v2 \
PREFLIGHT_ONLY=1 SLURM_JOB_ID=999603 \
bash slurm/run_model_training_603.sbatch \
  > "${STAGE}/validation/wrapper-preflight.stdout.txt" \
  2> "${STAGE}/validation/wrapper-preflight.stderr.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  config/campaigns/gpu_bundle_queue/061-model-training-603.json \
  slurm/run_model_training_603.sbatch \
  config/borzoi_grelu_converted_port_full_window_probe.json \
  scripts/probe_borzoi_grelu_converted_port_full_window.py \
  tests/unit/test_probe_borzoi_grelu_converted_port_full_window.py \
  scripts/gpu_bundle_dispatcher.py tests/unit/test_gpu_bundle_dispatcher.py \
  slurm/validate_model_training_603_queue_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
queue = json.loads((root / "validation/queue.json").read_text())
stderr = (root / "validation/model-tests.stderr.txt").read_text()
if (
    queue["enabled"] or not queue["ready"]
    or "Ran 4 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("model-training-603 queue admission receipt differs")
freeze_tree(root, {
    "artifact_class": "model_training_603_disabled_gpu_queue_admission",
    "bundle_id": "model-training-603", "family": "long_range_regulatory",
    "logical_tasks": 1, "queue_priority": 61, "queue_enabled": False,
    "queue_ready": True, "dispatcher_only": True, "qos": "nslab",
    "partition": "gpu", "gpu": "l40s:1", "cpus": 8,
    "memory_gb": 64, "wall_time": "08:00:00", "tests": 4,
    "biological_data_read": False, "development_outcomes_read": False,
    "sealed_sources_read": False, "training_performed": False,
    "prediction_performed": False, "gpu_job_submitted": False,
    "global_frozen_census_modified": False,
    "global_promotion_gate_modified": False,
    "queue_sha256": queue["queue_sha256"],
    "wrapper_sha256": queue["wrapper_sha256"], "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-training-603 disabled queue admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
