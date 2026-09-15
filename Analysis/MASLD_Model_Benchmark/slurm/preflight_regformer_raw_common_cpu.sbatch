#!/usr/bin/env bash
#SBATCH --job-name=model-check-031
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
ACQUISITION="${ROOT}/executions/regformer-acquisition-21066043"
INVENTORY="${ROOT}/executions/regformer-safe-inventory-21066069/checkpoint_inventory.json"
INSPECTION="${ROOT}/executions/regformer-weights-only-inspection-21066098/weights_only_inspection.json"
RUNTIME="${ROOT}/executions/environments/caduceus-native-torch2.2.0-21065466"
RUNTIME_SHA="3d59792d8a13ca10a2df3d6386bc266eb0e108def6279188793d1cfee4fc1167"
SOURCE_COMMIT="9618f1f4bb7b61c9b236d3ea8a6339c236f2da1b"
SCRIPT="${ROOT}/scripts/preflight_regformer_raw_common.py"
TEST="tests.unit.test_preflight_regformer_raw_common"
TARGET="${ROOT}/executions/model-check-031-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/model-check-031-${SLURM_JOB_ID}.XXXXXX")"
SOURCE="${WORK}/source"
RUNTIME_ENV="${WORK}/runtime"

preserve_failure() {
  status=$?
  if [[ ${status} -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed strict-restore evidence: %s\n' "${FAILED}" >&2
  fi
  if [[ -d "${WORK}" ]]; then
    rm -r "${WORK}"
  fi
  exit "${status}"
}
trap preserve_failure EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to replace preflight output: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}" "${RUNTIME_ENV}"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${ROOT}/src:${ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
CONTROL_PYTHON=python
"${CONTROL_PYTHON}" -m unittest "${TEST}" \
  > "${STAGE}/unit_tests.txt" 2>&1
printf '%s  %s\n' "${RUNTIME_SHA}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
"${CONTROL_PYTHON}" - "${RUNTIME}" <<'PY' \
  > "${STAGE}/runtime_artifact_verification.txt"
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("frozen runtime verified")
PY

git clone --no-checkout https://github.com/BGIResearch/RegFormer.git "${SOURCE}" \
  > "${STAGE}/git_clone.stdout.txt" 2> "${STAGE}/git_clone.stderr.txt"
git -C "${SOURCE}" checkout --detach "${SOURCE_COMMIT}" \
  > "${STAGE}/git_checkout.stdout.txt" 2> "${STAGE}/git_checkout.stderr.txt"
git -C "${SOURCE}" archive \
  --format=tar.gz \
  --prefix="RegFormer-${SOURCE_COMMIT}/" \
  --output="${STAGE}/regformer-source-${SOURCE_COMMIT}.tar.gz" \
  "${SOURCE_COMMIT}"
git -C "${SOURCE}" show --no-patch --format=fuller "${SOURCE_COMMIT}" \
  > "${STAGE}/source_commit.txt"

tar -xzf "${RUNTIME}/payload/caduceus-native-env.tar.gz" -C "${RUNTIME_ENV}"
"${RUNTIME_ENV}/bin/python" "${SCRIPT}" \
  --source "${SOURCE}" \
  --checkpoint "${ACQUISITION}/best_model.pt" \
  --vocab "${ACQUISITION}/vocab.json" \
  --graph "${ACQUISITION}/output_graph_no_cycles.dgl" \
  --config "${ACQUISITION}/pretrain.toml" \
  --inventory "${INVENTORY}" \
  --inspection "${INSPECTION}" \
  --runtime-artifacts "${RUNTIME}/ARTIFACTS.json" \
  --output "${STAGE}/strict_restore_receipt.json" \
  > "${STAGE}/strict_restore.stdout.txt" \
  2> "${STAGE}/strict_restore.stderr.txt"
"${RUNTIME_ENV}/bin/python" -m pip freeze --all \
  > "${STAGE}/runtime_pip_freeze.txt"
sha256sum \
  "${SCRIPT}" \
  "${ROOT}/tests/unit/test_preflight_regformer_raw_common.py" \
  "${ROOT}/slurm/preflight_regformer_raw_common_cpu.sbatch" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${ACQUISITION}/best_model.pt" \
  "${ACQUISITION}/vocab.json" \
  "${ACQUISITION}/output_graph_no_cycles.dgl" \
  "${ACQUISITION}/pretrain.toml" \
  "${INVENTORY}" \
  "${INSPECTION}" \
  "${STAGE}/regformer-source-${SOURCE_COMMIT}.tar.gz" \
  > "${STAGE}/source_and_input_sha256.txt"

"${CONTROL_PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "strict_restore_receipt.json").read_text())
if (
    receipt.get("status") != "pass_strict_restore_cpu_compatibility"
    or receipt.get("common_representation")
    != "raw_unrefined_avg_pooled_512d_encoder_embedding"
    or receipt.get("label_refined_embeddings_used") is not False
    or receipt.get("runtime", {}).get("strict_restore") is not True
    or receipt.get("runtime", {}).get("full_encoder_forward_executed") is not False
    or receipt.get("activation_disposition")
    != "blocked_pending_canonical_runtime_and_compiled_forward"
):
    raise SystemExit("RegFormer strict-restore receipt differs")
freeze_tree(
    root,
    {
        "model_id": "regformer",
        "artifact_class": "raw_common_cpu_strict_restore_preflight",
        "status": "pass_strict_restore_cpu_compatibility",
        "sealed_outcomes_read": False,
        "label_refined_embeddings_used": False,
        "activation_disposition": "blocked",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
rm -r "${WORK}"
printf 'RegFormer raw-common strict-restore candidate: %s\n' "${TARGET}"
