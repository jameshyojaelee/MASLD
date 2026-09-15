#!/usr/bin/env bash
#SBATCH --job-name=model-runtime-055
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE_ENV=/gpfs/commons/home/jameslee/micromamba/envs/scimilarity
AUDIT=${ROOT}/executions/scimilarity-existing-cellfm-runtime-audit-21066299
ARCHIVE=${ROOT}/executions/scimilarity-acquisition-21015744/model_v1.1.tar.gz
MODEL=${ROOT}/executions/scimilarity-model-20260823/model_v1.1
CONFIG=${ROOT}/config/artifacts/models/scimilarity
INVENTORY=${ROOT}/executions/scimilarity_v1_1_encoder-safe-inventory-21066203/checkpoint_inventory.json
WEIGHTS=${ROOT}/executions/scimilarity_v1_1_encoder-weights-only-inspection-21066254/weights_only_inspection.json
LICENSE=${SOURCE_ENV}/lib/python3.10/site-packages/scimilarity-0.4.1.dist-info/licenses/LICENSE
NN_MODELS=${SOURCE_ENV}/lib/python3.10/site-packages/scimilarity/nn_models.py
TARGET=${ROOT}/executions/scimilarity-v1-1-compat-runtime-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
WORK=$(mktemp -d "${TMPDIR:-/tmp}/scimilarity-runtime.${SLURM_JOB_ID}.XXXXXX")

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
  printf 'Refusing to overwrite SCimilarity compatibility runtime.\n' >&2
  exit 1
fi
if [[ ! -x "${SOURCE_ENV}/bin/python" ]]; then
  printf 'SCimilarity source runtime is absent.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/payload" "${STAGE}/validation" "${WORK}/runtime"
cd "${ROOT}"
printf '%s  %s\n' \
  2437806f5ab0c32a9ba74760545897a47de0a701133ea840d2209a835a64926f "${AUDIT}/runtime_audit.json" \
  afd3db4500c816f3a6698bf7af67881d1fc45707b661e041e29c1251781e5939 "${AUDIT}/pip_freeze.txt" \
  1d7e8b2e174da474eedd352b99ac2fb8ea544e1f26a8b4ed6bad9de57a08f0e8 "${AUDIT}/conda_explicit.txt" \
  3e1a15188a0092cf0924a82febd6a44e373975712541b66b7a7f98d0859557cd "${CONFIG}/checkpoints.json" \
  c43c0df5c2d39703ad09ff81aebafc8c21f7a036d4738a2e605c0ff85d983b01 "${CONFIG}/exposure_audit.json" \
  c0e83d7298c25025cefd313eb05e68979d6e8a77d03e3bad04d29118cdab9abf "${CONFIG}/development_crosswalk.json" \
  ec597b0b7b193e17d5ff2db891afd05cc4cc41dc53ce384d6712dc32646eca36 "${INVENTORY}" \
  26720d7dfa0276617bc744fb390362321750246cfe53827bf3db51053b51eeeb "${WEIGHTS}" \
  a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f "${MODEL}/encoder.ckpt" \
  0a9aeacf7fc19399426273d8cfab8683253ec56b51be556d42edd115f125b8ca "${MODEL}/gene_order.tsv" \
  938fb579905acc3f45d10f20cd18da74f2e754300220f9669d0e5c5698f85cac "${LICENSE}" \
  5b3edf1a57ae653202d9655e049b516ec31cfa14e81daaab7c96b647fec6c110 "${NN_MODELS}" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
CONTROL_PYTHON=$(command -v python)
PYTHONPATH="${ROOT}:${ROOT}/src" "${CONTROL_PYTHON}" scripts/audit_scimilarity_v1_1_activation.py \
  --archive "${ARCHIVE}" \
  --model-root "${MODEL}" \
  --code-license "${LICENSE}" \
  --nn-models "${NN_MODELS}" \
  --checkpoint-config "${CONFIG}/checkpoints.json" \
  --exposure-audit "${CONFIG}/exposure_audit.json" \
  --crosswalk "${CONFIG}/development_crosswalk.json" \
  --output "${STAGE}/activation_audit.json" \
  > "${STAGE}/validation/activation-audit.stdout.json"

export PYTHONPATH="${ROOT}:${ROOT}/src"
"${SOURCE_ENV}/bin/python" - "${MODEL}/encoder.ckpt" "${STAGE}/validation/source_embedding.npy" <<'PY'
from pathlib import Path
import sys
import numpy as np
from scipy import sparse

from scripts.extract_scimilarity_v1_1_frozen_screen_50000 import _load_encoder, embed_batch

model, metadata = _load_encoder(Path(sys.argv[1]), device="cpu")
matrix = sparse.csr_matrix((np.asarray([7.5, 2.5], dtype=np.float32), ([0, 1], [0, 1])), shape=(2, 28231))
first = embed_batch(model, matrix, device="cpu")
second = embed_batch(model, matrix, device="cpu")
if not np.array_equal(first, second):
    raise SystemExit("source-runtime deterministic repeat differs")
np.save(sys.argv[2], first, allow_pickle=False)
PY

tar --sort=name --mtime='UTC 2026-08-24' --owner=0 --group=0 --numeric-owner \
  --pax-option=delete=atime,delete=ctime -C "${SOURCE_ENV}" -cf - . \
  | gzip -1n > "${STAGE}/payload/scimilarity-v1-1-compat-env.tar.gz"
sha256sum "${STAGE}/payload/scimilarity-v1-1-compat-env.tar.gz" \
  > "${STAGE}/validation/runtime-archive.sha256"
tar -xzf "${STAGE}/payload/scimilarity-v1-1-compat-env.tar.gz" -C "${WORK}/runtime"
"${WORK}/runtime/bin/python" -m unittest -v tests.unit.test_scimilarity_v1_1_frozen_screen_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${WORK}/runtime/bin/python" - "${MODEL}/encoder.ckpt" \
  "${STAGE}/validation/source_embedding.npy" "${STAGE}/runtime_receipt.json" <<'PY'
from pathlib import Path
import json
import sys
import numpy as np
from scipy import sparse
import torch

from scripts.extract_scimilarity_v1_1_frozen_screen_50000 import _load_encoder, embed_batch

checkpoint = Path(sys.argv[1])
expected = np.load(sys.argv[2], allow_pickle=False)
model, loading = _load_encoder(checkpoint, device="cpu")
matrix = sparse.csr_matrix((np.asarray([7.5, 2.5], dtype=np.float32), ([0, 1], [0, 1])), shape=(2, 28231))
observed = embed_batch(model, matrix, device="cpu")
if not np.array_equal(expected, observed):
    raise SystemExit("relocated runtime numeric parity differs")
receipt = {
    "schema_version": "masld-bench-scimilarity-v1-1-compat-runtime-v1",
    "status": "pass_cpu_runtime_strict_restore_and_numeric_parity",
    "model_id": "scimilarity_v1_1",
    "checkpoint_sha256": "a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f",
    "strict_weights_only_restore": True,
    "architecture_source_sha256": "5b3edf1a57ae653202d9655e049b516ec31cfa14e81daaab7c96b647fec6c110",
    "parameter_count": loading["parameter_count"],
    "embedding_width": 128,
    "numeric_parity_max_abs_diff": float(np.max(np.abs(expected - observed))),
    "torch_version": torch.__version__,
    "raw_unrefined_embeddings_only": True,
    "released_reference_index_used": False,
    "downstream_head_fit": False,
    "evaluation_label_columns_read": [],
    "sealed_outcomes_read": False,
}
Path(sys.argv[3]).write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
PY
"${WORK}/runtime/bin/python" -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
cp "${AUDIT}/runtime_audit.json" "${STAGE}/validation/source-runtime-audit.json"
cp "${AUDIT}/conda_explicit.txt" "${STAGE}/validation/source-conda-explicit.txt"
cp "${AUDIT}/pip_freeze.txt" "${STAGE}/validation/source-pip-freeze.txt"
sha256sum \
  scripts/audit_scimilarity_v1_1_activation.py \
  scripts/build_scimilarity_v1_1_frozen_screen_50000.py \
  scripts/extract_scimilarity_v1_1_frozen_screen_50000.py \
  tests/unit/test_scimilarity_v1_1_frozen_screen_50000.py \
  slurm/build_scimilarity_v1_1_compat_runtime_cpu.sbatch \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
runtime = json.loads((root / "runtime_receipt.json").read_text(encoding="utf-8"))
audit = json.loads((root / "activation_audit.json").read_text(encoding="utf-8"))
if (
    runtime.get("status") != "pass_cpu_runtime_strict_restore_and_numeric_parity"
    or runtime.get("checkpoint_sha256") != "a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f"
    or runtime.get("embedding_width") != 128
    or runtime.get("numeric_parity_max_abs_diff") != 0.0
    or runtime.get("strict_weights_only_restore") is not True
    or runtime.get("released_reference_index_used") is not False
    or runtime.get("downstream_head_fit") is not False
    or runtime.get("evaluation_label_columns_read")
    or runtime.get("sealed_outcomes_read") is not False
    or audit.get("status") != "pass_exact_assets_terms_and_exposure"
    or audit.get("central_registry_reconciliation_required_before_campaign_admission") is not True
):
    raise SystemExit("SCimilarity compatibility runtime receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "scimilarity_v1_1_compatibility_runtime",
        "model_id": "scimilarity_v1_1",
        "checkpoint_sha256": runtime["checkpoint_sha256"],
        "runtime_archive_sha256": sha256_file(root / "payload/scimilarity-v1-1-compat-env.tar.gz"),
        "runtime_claim": "compatibility_snapshot_not_upstream_canonical_environment",
        "strict_weights_only_restore": True,
        "cpu_numeric_parity": True,
        "released_reference_index_used": False,
        "downstream_head_fit": False,
        "evaluation_labels_read": False,
        "sealed_outcomes_read": False,
        "status": "passed_cpu_gate",
    },
)
PY
trap - EXIT
rm -rf -- "${WORK}"
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
