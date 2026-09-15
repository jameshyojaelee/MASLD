#!/usr/bin/env bash
#SBATCH --job-name=model-check-056
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
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/scimilarity-v1-1-compat-runtime-21079070
RUNTIME_SHA256=8f66e2e478742665072a23ed4b115346830fcee516094c59f20164c94b3a5a4a
RUNTIME_ARCHIVE_SHA256=e2d0a97c8631d3cdacb4234fa7bc167ddd8fc391c3dd2a2b7dd23ebaebe88548
FIXTURE=${ROOT}/executions/scimilarity-v1-1-frozen-screen-fixture-21079292
FIXTURE_SHA256=b6d50557107b5e6a878480741ec6a4c781663218876943ae205578b1a4ecd1e5
PREFLIGHT=${ROOT}/executions/scimilarity-v1-1-frozen-screen-preflight-21079643
PREFLIGHT_SHA256=8dcf5c8d91f59ebfa5cbb9d55be7d3299c0eeaae6522b9be7c3a893306eeb6e9
STAGED_QUEUE=${ROOT}/config/campaigns/gpu_bundle_queue_staged/054-model-training-504.json
STAGED_QUEUE_SHA256=365e35cfa46bfb12509f93ea5ba3d4f1d142a63c1cd54dfabe2cd4e78307fbd4
LIVE_QUEUE=${ROOT}/config/campaigns/gpu_bundle_queue/054-model-training-504.json
WRAPPER=${ROOT}/slurm/run_scimilarity_v1_1_frozen_screen_50000_resumable.sbatch
WRAPPER_SHA256=e35a7f5b5c6d7bbce6ffea53bda3fcd237db87658ee5d43bbd01fbb9938fc878
EXTRACTOR=${ROOT}/scripts/extract_scimilarity_v1_1_frozen_screen_50000_resumable.py
EXTRACTOR_SHA256=e664cef093641202dd43989bd284e3c126ec644262913a39f6bc5d67ad0b15bd
TEST=${ROOT}/tests/unit/test_scimilarity_v1_1_resumable_extraction.py
TEST_SHA256=a0adbaa75b5d5e0e6edafd610f784620ac7ba91b3c79b609fd397e3b0bbde44c
TARGET=${ROOT}/executions/model-training-504-queue-validation-v1
STAGE=${TARGET}.staging-${SLURM_JOB_ID}
FAILED=${TARGET}.failed-${SLURM_JOB_ID}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/model-check-056.${SLURM_JOB_ID}.XXXXXX")

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  if [[ -d "${WORK}" ]]; then
    rm -rf -- "${WORK}"
  fi
  exit "${exit_code}"
}
trap finish EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite SCimilarity GPU queue validation.\n' >&2
  exit 1
fi
if [[ -e "${LIVE_QUEUE}" ]]; then
  printf 'Live priority-54 queue item exists before validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${WORK}/runtime"
cd "${ROOT}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${RUNTIME_ARCHIVE_SHA256}" "${RUNTIME}/payload/scimilarity-v1-1-compat-env.tar.gz" \
  "${FIXTURE_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  "${PREFLIGHT_SHA256}" "${PREFLIGHT}/ARTIFACTS.json" \
  "${STAGED_QUEUE_SHA256}" "${STAGED_QUEUE}" \
  "${WRAPPER_SHA256}" "${WRAPPER}" \
  "${EXTRACTOR_SHA256}" "${EXTRACTOR}" \
  "${TEST_SHA256}" "${TEST}" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"
bash -n "${WRAPPER}"

module purge
module load python/3.11.5-GCCcore-13.2.0
CONTROL_PYTHON=$(command -v python)
PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - \
  "${RUNTIME}" "${FIXTURE}" "${PREFLIGHT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
tar -xzf "${RUNTIME}/payload/scimilarity-v1-1-compat-env.tar.gz" -C "${WORK}/runtime"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${WORK}/runtime/bin/python" -m unittest -v \
  tests.unit.test_gpu_bundle_dispatcher \
  tests.unit.test_scimilarity_v1_1_frozen_screen_50000 \
  tests.unit.test_scimilarity_v1_1_resumable_extraction \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${WORK}/runtime/bin/python" - "${STAGED_QUEUE}" "${ROOT}" "${WRAPPER}" \
  > "${STAGE}/validation/prepublication.json" <<'PY'
from pathlib import Path
import json
import sys

from scripts.gpu_bundle_dispatcher import load_item

queue = Path(sys.argv[1])
root = Path(sys.argv[2]).resolve(strict=True)
wrapper = Path(sys.argv[3]).read_text(encoding="utf-8")
item = load_item(queue, root)
required_header = {
    "#SBATCH --job-name=model-training-504",
    "#SBATCH --partition=gpu",
    "#SBATCH --account=nslab",
    "#SBATCH --qos=nslab",
    "#SBATCH --gres=gpu:l40s:1",
    "#SBATCH --cpus-per-task=4",
    "#SBATCH --mem=64G",
    "#SBATCH --time=04:00:00",
    "#SBATCH --requeue",
    "#SBATCH --signal=B:USR1@120",
}
if (
    item["_ready"]
    or item["bundle_id"] != "model-training-504"
    or item["priority"] != 54
    or item["logical_tasks"] != 1
    or item["family"] != "cell_foundation"
    or not required_header <= set(wrapper.splitlines())
    or "#SBATCH --array" in wrapper
    or "--qos=innovation" in wrapper
    or "fit_predict_common_cell_heads" in wrapper
    or "evaluate_" in wrapper
):
    raise SystemExit("staged SCimilarity queue contract differs")
print(
    json.dumps(
        {
            "bundle_id": item["bundle_id"],
            "priority": item["priority"],
            "ready_before_validation": False,
            "live_queue_mutated": False,
        },
        sort_keys=True,
    )
)
PY
"${WORK}/runtime/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  scripts/gpu_bundle_dispatcher.py \
  scripts/extract_scimilarity_v1_1_frozen_screen_50000.py \
  scripts/extract_scimilarity_v1_1_frozen_screen_50000_resumable.py \
  tests/unit/test_gpu_bundle_dispatcher.py \
  tests/unit/test_scimilarity_v1_1_frozen_screen_50000.py \
  tests/unit/test_scimilarity_v1_1_resumable_extraction.py \
  slurm/run_scimilarity_v1_1_frozen_screen_50000_resumable.sbatch \
  slurm/validate_scimilarity_v1_1_gpu_queue_cpu.sbatch \
  config/campaigns/gpu_bundle_queue_staged/054-model-training-504.json \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="${ROOT}:${ROOT}/src" "${CONTROL_PYTHON}" - "${STAGE}" \
  "${STAGED_QUEUE_SHA256}" "${WRAPPER_SHA256}" \
  "${EXTRACTOR_SHA256}" "${TEST_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = {
    "schema_version": "masld-bench-scimilarity-v1-1-gpu-queue-validation-v1",
    "status": "passed",
    "bundle_id": "model-training-504",
    "priority": 54,
    "logical_tasks": 1,
    "staged_queue_sha256": sys.argv[2],
    "wrapper_sha256": sys.argv[3],
    "resumable_extractor_sha256": sys.argv[4],
    "resumability_test_sha256": sys.argv[5],
    "partition": "gpu",
    "account": "nslab",
    "qos": "nslab",
    "gres": "gpu:l40s:1",
    "cpus": 4,
    "memory": "64G",
    "time": "04:00:00",
    "array": False,
    "job_name": "model-training-504",
    "dispatcher_only": True,
    "shared_gate": {"cap": 5, "max_running": 4, "max_pending": 1},
    "atomic_chunk_resume_tested": True,
    "atomic_final_publication": True,
    "raw_unrefined_embeddings_only": True,
    "released_reference_index_used": False,
    "metrics_run": False,
    "head_fit": False,
    "evaluation_labels_read": False,
    "sealed_outcomes_read": False,
    "gpu_submitted": False,
    "live_queue_mutated": False,
}
(stage / "validation/receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
freeze_tree(
    stage,
    {
        "artifact_class": "scimilarity_v1_1_gpu_queue_contract_validation",
        "bundle_id": "model-training-504",
        "priority": 54,
        "staged_queue_sha256": sys.argv[2],
        "wrapper_sha256": sys.argv[3],
        "gpu_submitted": False,
        "live_queue_mutated": False,
        "status": "passed",
    },
)
PY
PYTHONPATH="${ROOT}:${ROOT}/src" "${CONTROL_PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import publish_directory_noreplace

publish_directory_noreplace(Path(sys.argv[1]), Path(sys.argv[2]))
PY
PYTHONPATH="${ROOT}:${ROOT}/src" "${CONTROL_PYTHON}" - "${STAGED_QUEUE}" "${ROOT}" "${TARGET}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree
from scripts.gpu_bundle_dispatcher import load_item

verify_frozen_tree(Path(sys.argv[3]))
item = load_item(Path(sys.argv[1]), Path(sys.argv[2]).resolve(strict=True))
if not item["_ready"] or item["bundle_id"] != "model-training-504":
    raise SystemExit("staged SCimilarity queue is not ready after validation")
print(item["_queue_sha256"])
PY
if [[ -e "${LIVE_QUEUE}" ]]; then
  printf 'Live priority-54 queue item changed during validation.\n' >&2
  exit 1
fi

trap - EXIT
rm -rf -- "${WORK}"
sha256sum "${TARGET}/ARTIFACTS.json"
