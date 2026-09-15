#!/usr/bin/env bash
#SBATCH --job-name=model-check-265
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
QUEUE="${ROOT}/config/campaigns/gpu_bundle_queue/060-model-training-602.json"
TARGET="${ROOT}/executions/model-check-265-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-training-602 queue validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "7f28ec8d5db5bba6ff755ed2d39fc26dbfd90026b629ed2a646c48bcf4b07e9e ${QUEUE}" \
  "7cc5201dca83652fc56dab912514ec9ca6d85f225d04ce65a8e4ff29dd6ed04b ${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "06eb07080d4cae46eeb373fd309dadc971e149a74ce857adea1c3945543e901e ${ROOT}/tests/unit/test_gpu_bundle_dispatcher.py" \
  "af7b86cc633700458eaf68ec8e567e30a467f043bbafbad4482d5300341f542f ${ROOT}/slurm/run_model_training_602.sbatch" \
  "951f1349043dd9f8c2c5697bd7f64d93fa6f16b1f28b6dc337b1ef420d82f904 ${ROOT}/scripts/probe_corgi_regular_film_head_l40s.py" \
  "2259363a1f8d0729f095c904f1f59ced21c349aa8ec140b98a1445f30cf91dfd ${ROOT}/tests/unit/test_probe_corgi_regular_film_head_l40s.py" \
  "6dd6242696f458efbec72bf73a36938802487a6d0e8f4e1ae2a401a3872129ca ${ROOT}/src/masld_bench/adapters/corgi_film_head.py" \
  "368585b968f690646ccd60a59df146c24aaf429b28058225491c1c991d3fe596 ${ROOT}/tests/unit/test_corgi_film_head_adapter.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_gpu_bundle_dispatcher \
  > "${STAGE}/validation/dispatcher-tests.stdout.txt" \
  2> "${STAGE}/validation/dispatcher-tests.stderr.txt"
PYTHONPATH="${ROOT}" python - "${QUEUE}" "${ROOT}" \
  > "${STAGE}/validation/queue.json" <<'PY'
import json
from pathlib import Path
import sys
from scripts.gpu_bundle_dispatcher import load_item

item = load_item(Path(sys.argv[1]), Path(sys.argv[2]))
if (
    item["bundle_id"] != "model-training-602"
    or item["priority"] != 60
    or item["enabled"]
    or not item["_ready"]
    or item["logical_tasks"] != 1
    or item["family"] != "context_regulatory"
):
    raise SystemExit("disabled model-training-602 queue item differs")
print(json.dumps({
    "bundle_id": item["bundle_id"], "priority": item["priority"],
    "enabled": item["enabled"], "ready": item["_ready"],
    "queue_sha256": item["_queue_sha256"],
    "wrapper_sha256": item["wrapper_sha256"],
}, indent=2, sort_keys=True))
PY
MASLD_GPU_DISPATCHER_AUTHORITY=sequence_gpu_dispatcher_v2 \
PREFLIGHT_ONLY=1 SLURM_JOB_ID=999602 \
bash slurm/run_model_training_602.sbatch \
  > "${STAGE}/validation/wrapper-preflight.stdout.txt" \
  2> "${STAGE}/validation/wrapper-preflight.stderr.txt"
module purge
module load singularity/4.2.2
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --pwd "${ROOT}" --env "PYTHONPATH=${ROOT}" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v \
  tests.unit.test_corgi_film_head_adapter \
  tests.unit.test_probe_corgi_regular_film_head_l40s \
  > "${STAGE}/validation/model-tests.stdout.txt" \
  2> "${STAGE}/validation/model-tests.stderr.txt"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --pwd "${ROOT}" "${IMAGE}" \
  python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
cd "${ROOT}"
sha256sum config/campaigns/gpu_bundle_queue/060-model-training-602.json \
  scripts/gpu_bundle_dispatcher.py tests/unit/test_gpu_bundle_dispatcher.py \
  slurm/run_model_training_602.sbatch scripts/probe_corgi_regular_film_head_l40s.py \
  tests/unit/test_probe_corgi_regular_film_head_l40s.py \
  src/masld_bench/adapters/corgi_film_head.py \
  tests/unit/test_corgi_film_head_adapter.py \
  slurm/validate_model_training_602_queue_cpu.sbatch > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
queue = json.loads((root / "validation/queue.json").read_text())
model_stderr = (root / "validation/model-tests.stderr.txt").read_text()
if (
    queue["enabled"] or not queue["ready"]
    or "Ran 8 tests" not in model_stderr
    or not model_stderr.rstrip().endswith("OK")
):
    raise SystemExit("model-training-602 queue admission receipt differs")
freeze_tree(root, {
    "artifact_class": "model_training_602_disabled_gpu_queue_admission",
    "bundle_id": "model-training-602", "family": "context_regulatory",
    "logical_tasks": 1, "queue_priority": 60, "queue_enabled": False,
    "queue_ready": True, "dispatcher_only": True, "qos": "nslab",
    "partition": "gpu", "gpu": "l40s:1", "cpus": 8,
    "memory_gb": 64, "wall_time": "08:00:00", "tests": 8,
    "development_atac_outcomes_read": False,
    "histology_or_disease_labels_read": False,
    "test_or_sealed_features_or_outcomes_read": False,
    "training_performed": False, "prediction_performed": False,
    "gpu_job_submitted": False, "global_census_modified": False,
    "global_promotion_gate_modified": False,
    "queue_sha256": queue["queue_sha256"],
    "wrapper_sha256": queue["wrapper_sha256"], "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-training-602 disabled queue admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
