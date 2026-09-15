#!/usr/bin/env bash
#SBATCH --job-name=model-check-076
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
EMBEDDINGS=${ROOT}/executions/model-training-500-21080428
EMBEDDINGS_SHA256=fc4601988491b086d5e89f0ff35f2f73d9c27cf74846a6ff7895efc130c825a5
QUEUE=${ROOT}/config/campaigns/gpu_bundle_queue/052-model-training-505.json
QUEUE_SHA256=64e2ad1c8985b22dea6299e24ac9abe163080fa66d91632f1b050a52e36e0e22
WRAPPER=${ROOT}/slurm/run_common_cell_head_bundle_50000.sbatch
WRAPPER_SHA256=40af0f57210120ac50c1ca19f138a44f865aa03346af02dc56462f4e2b1cb3cf
TARGET=${ROOT}/executions/transcriptformer-common-head-queue-validation-v1
STAGE=${TARGET}.staging-${SLURM_JOB_ID}
FAILED=${TARGET}.failed-${SLURM_JOB_ID}

preserve_failure() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${exit_code}"
}
trap preserve_failure EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer common-head queue validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  "${EMBEDDINGS_SHA256}" "${EMBEDDINGS}/ARTIFACTS.json" \
  "${QUEUE_SHA256}" "${QUEUE}" \
  "${WRAPPER_SHA256}" "${WRAPPER}" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"
bash -n "${WRAPPER}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_gpu_bundle_dispatcher \
  tests.unit.test_fit_predict_common_cell_heads_study_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${RUNTIME}/env/bin/python" - "${QUEUE}" "${ROOT}" \
  > "${STAGE}/validation/prepublication.json" <<'PY'
from pathlib import Path
import json
import sys

from scripts.gpu_bundle_dispatcher import load_item

item = load_item(Path(sys.argv[1]), Path(sys.argv[2]).resolve(strict=True))
if (
    item["_ready"]
    or item["bundle_id"] != "model-training-505"
    or item["priority"] != 52
    or item["logical_tasks"] != 30
    or item["exports"].get("MASLD_GPU_DISPATCHER_AUTHORITY")
    != "sequence_gpu_dispatcher_v2"
):
    raise SystemExit("TranscriptFormer common-head queue prepublication state differs")
print(json.dumps({"bundle_id": item["bundle_id"], "ready_before_validation": False}, sort_keys=True))
PY
"${RUNTIME}/env/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/gpu_bundle_dispatcher.py \
  scripts/fit_predict_common_cell_heads_study_50000.py \
  tests/unit/test_gpu_bundle_dispatcher.py \
  tests/unit/test_fit_predict_common_cell_heads_study_50000.py \
  slurm/run_common_cell_head_bundle_50000.sbatch \
  slurm/validate_transcriptformer_common_head_gpu_queue.sbatch \
  config/campaigns/gpu_bundle_queue/052-model-training-505.json \
  > "${STAGE}/validation/source.sha256"
"${RUNTIME}/env/bin/python" - "${STAGE}" "${QUEUE_SHA256}" "${WRAPPER_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
tests = (stage / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if "OK" not in tests:
    raise SystemExit("TranscriptFormer common-head queue tests did not pass")
freeze_tree(
    stage,
    {
        "artifact_class": "transcriptformer_common_head_gpu_queue_validation",
        "bundle_id": "model-training-505",
        "queue_sha256": sys.argv[2],
        "wrapper_sha256": sys.argv[3],
        "logical_tasks": 30,
        "dispatcher_only": True,
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
