#!/usr/bin/env bash
#SBATCH --job-name=model-test-208
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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
QUEUE=${ROOT}/config/campaigns/gpu_bundle_queue/053-model-work-208.json
QUEUE_SHA256=3bdd62df02d821b1e4355876a8ba2544dc69baf656806648b6309b4c7e783ff7
WRAPPER=${ROOT}/slurm/run_uce_frozen_screen_50000.sbatch
WRAPPER_SHA256=09958f9c691144dbfe92723dc9c9a75a9c014e8c36e7dd868027f3f635921c94
FIXTURE=${ROOT}/executions/model-preflight-507-21079352
FIXTURE_SHA256=6b9e513ef171dbe9b22ede423fe008c333dbb271458e1a662f7c6684787c3260
TARGET=${ROOT}/executions/model-work-208-queue-validation-v1
STAGE=${TARGET}.staging-${SLURM_JOB_ID}
FAILED=${TARGET}.failed-${SLURM_JOB_ID}

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${exit_code}"
}
trap finish EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite UCE queue validation output.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${QUEUE_SHA256}" "${QUEUE}" \
  "${WRAPPER_SHA256}" "${WRAPPER}" \
  "${FIXTURE_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"
bash -n "${WRAPPER}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_gpu_bundle_dispatcher \
  tests.unit.test_uce_gpu_queue_contract \
  tests.unit.test_uce_frozen_screen_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${RUNTIME}/env/bin/python" - "${QUEUE}" "${ROOT}" \
  > "${STAGE}/validation/prepublication.json" <<'PY'
from pathlib import Path
import json
import sys
from scripts.gpu_bundle_dispatcher import load_item

item = load_item(Path(sys.argv[1]), Path(sys.argv[2]).resolve(strict=True))
if item["_ready"]:
    raise SystemExit("UCE queue became ready before validation publication")
print(json.dumps({"bundle_id": item["bundle_id"], "ready_before_validation": False}))
PY
sha256sum \
  scripts/gpu_bundle_dispatcher.py \
  scripts/uce_frozen_screen_50000_common.py \
  scripts/build_uce_frozen_screen_50000.py \
  scripts/extract_uce_frozen_screen_50000.py \
  slurm/run_uce_frozen_screen_50000.sbatch \
  slurm/test_uce_gpu_queue_contract_cpu.sbatch \
  tests/unit/test_gpu_bundle_dispatcher.py \
  tests/unit/test_uce_gpu_queue_contract.py \
  tests/unit/test_uce_frozen_screen_50000.py \
  config/campaigns/gpu_bundle_queue/053-model-work-208.json \
  config/artifacts/models/uce/frozen_screen_50000.json \
  > "${STAGE}/validation/source.sha256"
"${RUNTIME}/env/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/pip-freeze.txt"
"${RUNTIME}/env/bin/python" - "${STAGE}" "${QUEUE_SHA256}" \
  "${WRAPPER_SHA256}" "${FIXTURE_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = {
    "schema_version": "masld-bench-uce-gpu-queue-validation-v1",
    "status": "passed",
    "bundle_id": "model-work-208",
    "priority": 53,
    "logical_tasks": ["uce_4l", "uce_33l"],
    "execution_order": ["uce_4l", "uce_33l"],
    "queue_sha256": sys.argv[2],
    "wrapper_sha256": sys.argv[3],
    "fixture_artifacts_sha256": sys.argv[4],
    "partition": "gpu",
    "account": "nslab",
    "qos": "nslab",
    "gres": "gpu:l40s:1",
    "cpus": 12,
    "memory": "128G",
    "time": "120:00:00",
    "array": False,
    "job_name": "model-work-208",
    "dispatcher_only": True,
    "shared_gate": {"cap": 5, "max_running": 4, "max_pending": 1},
    "completed_model_chunk_safe_skip": True,
    "separate_immutable_output_roots": True,
    "metrics_run": False,
    "head_fit": False,
    "manual_gpu_submission": False,
}
(stage / "validation/receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
freeze_tree(
    stage,
    {
        "artifact_class": "uce_gpu_queue_contract_validation",
        "bundle_id": "model-work-208",
        "priority": 53,
        "queue_sha256": sys.argv[2],
        "wrapper_sha256": sys.argv[3],
        "gpu_submitted": False,
        "status": "passed",
    },
)
PY
"${RUNTIME}/env/bin/python" - "${STAGE}" "${TARGET}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import publish_directory_noreplace

publish_directory_noreplace(Path(sys.argv[1]), Path(sys.argv[2]))
PY
"${RUNTIME}/env/bin/python" - "${QUEUE}" "${ROOT}" "${TARGET}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
from scripts.gpu_bundle_dispatcher import load_item

verify_frozen_tree(Path(sys.argv[3]))
item = load_item(Path(sys.argv[1]), Path(sys.argv[2]).resolve(strict=True))
if not item["_ready"] or item["bundle_id"] != "model-work-208":
    raise SystemExit("UCE queue is not ready after immutable validation publication")
print(item["_queue_sha256"])
PY
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
