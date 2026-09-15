#!/usr/bin/env bash
#SBATCH --job-name=model-validate-207
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/geneformer-v2-316m-common-head-queue-validation-v1
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Geneformer common-head queue validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_prepare_common_cell_embeddings_50000 \
  tests.unit.test_audit_evaluate_common_cell_heads_study_50000 \
  tests.unit.test_fit_predict_common_cell_heads_study_50000 \
  tests.unit.test_geneformer_v2_316m_common_head_queue \
  tests.unit.test_gpu_bundle_dispatcher \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" - <<'PY' > "${STAGE}/validation/admission.json"
from pathlib import Path
import json

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.hashing import sha256_file
from scripts.gpu_bundle_dispatcher import load_item

root = Path.cwd().resolve(strict=True)
queue = root / "config/campaigns/gpu_bundle_queue/058-model-training-506.json"
item = load_item(queue, root)
missing = [
    path
    for path in item["required_paths"]
    if not (root / path).is_file()
]
gate = "executions/geneformer-v2-316m-common-head-queue-validation-v1/ARTIFACTS.json"
if (
    item["bundle_id"] != "model-training-506"
    or item["logical_tasks"] != 30
    or item["exports"]["MASLD_GPU_DISPATCHER_AUTHORITY"]
    != "sequence_gpu_dispatcher_v2"
    or item["exports"]["EMBEDDINGS_ARTIFACTS_SHA256"]
    != "a675daa170020f4e12eac2eb6dcf046cac4304511789805f802e56d08900f86d"
    or missing != [gate]
):
    raise SystemExit("Geneformer common-head dispatcher admission differs")
verify_frozen_tree(root / "executions/model-data-053-21088295/embeddings")
print(
    json.dumps(
        {
            "bundle_id": item["bundle_id"],
            "logical_tasks": item["logical_tasks"],
            "queue_sha256": item["_queue_sha256"],
            "wrapper_sha256": item["wrapper_sha256"],
            "only_missing_path_is_self_gate": missing == [gate],
            "status": "pass_fail_closed_until_validation_freeze",
        },
        sort_keys=True,
    )
)
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/prepare_common_cell_embeddings_50000.py \
  scripts/fit_predict_common_cell_heads_study_50000.py \
  scripts/audit_aggregate_common_cell_heads_study_50000.py \
  scripts/score_common_cell_heads_study_50000.py \
  scripts/gpu_bundle_dispatcher.py \
  tests/unit/test_prepare_common_cell_embeddings_50000.py \
  tests/unit/test_fit_predict_common_cell_heads_study_50000.py \
  tests/unit/test_audit_evaluate_common_cell_heads_study_50000.py \
  tests/unit/test_geneformer_v2_316m_common_head_queue.py \
  tests/unit/test_gpu_bundle_dispatcher.py \
  slurm/prepare_geneformer_v2_316m_common_50000_cpu.sbatch \
  slurm/run_common_cell_head_bundle_50000.sbatch \
  slurm/audit_aggregate_common_cell_heads_study_50000.sbatch \
  slurm/score_common_cell_heads_study_50000.sbatch \
  slurm/validate_geneformer_v2_316m_common_head_queue_cpu.sbatch \
  config/campaigns/gpu_bundle_queue/058-model-training-506.json \
  config/campaigns/geneformer_v2_316m_common_head_chain_v1.json \
  executions/model-data-053-21088295/ARTIFACTS.json \
  executions/model-data-053-21088295/embeddings/ARTIFACTS.json \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
admission = json.loads((root / "validation/admission.json").read_text())
if (
    admission["status"] != "pass_fail_closed_until_validation_freeze"
    or not admission["only_missing_path_is_self_gate"]
    or admission["logical_tasks"] != 30
):
    raise SystemExit("Geneformer common-head admission receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "geneformer_v2_316m_common_head_queue_validation",
        "model_id": "geneformer_v2_316m",
        "bundle_id": "model-training-506",
        "logical_tasks": 30,
        "dispatcher_only": True,
        "gpu_account": "nslab",
        "gpu_qos": "nslab",
        "gpu_array": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
