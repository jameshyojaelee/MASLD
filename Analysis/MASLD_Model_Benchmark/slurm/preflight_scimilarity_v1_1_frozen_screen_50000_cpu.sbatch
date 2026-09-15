#!/usr/bin/env bash
#SBATCH --job-name=model-check-055
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/scimilarity-v1-1-compat-runtime-21079070
RUNTIME_SHA256=8f66e2e478742665072a23ed4b115346830fcee516094c59f20164c94b3a5a4a
RUNTIME_ARCHIVE_SHA256=e2d0a97c8631d3cdacb4234fa7bc167ddd8fc391c3dd2a2b7dd23ebaebe88548
FIXTURE=${ROOT}/executions/scimilarity-v1-1-frozen-screen-fixture-21079292
FIXTURE_SHA256=b6d50557107b5e6a878480741ec6a4c781663218876943ae205578b1a4ecd1e5
ACQUISITION=${ROOT}/executions/scimilarity-acquisition-21015744
ACQUISITION_RECORD_SHA256=f0643d5a909b2e8f6b5aaebe47af95ced985cb17a50829cd0e76147f238acebb
MODEL=${ROOT}/executions/scimilarity-model-20260823/model_v1.1
CONFIG=${ROOT}/config/artifacts/models/scimilarity
MODEL_REGISTRY=${ROOT}/config/models/cell_foundation.toml
TARGET=${ROOT}/executions/scimilarity-v1-1-frozen-screen-preflight-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
WORK=$(mktemp -d "${TMPDIR:-/tmp}/scimilarity-50k-preflight.${SLURM_JOB_ID}.XXXXXX")

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
  printf 'Refusing to overwrite SCimilarity 50k CPU preflight.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${WORK}/runtime"
cd "${ROOT}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${RUNTIME_ARCHIVE_SHA256}" "${RUNTIME}/payload/scimilarity-v1-1-compat-env.tar.gz" \
  "${FIXTURE_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  "${ACQUISITION_RECORD_SHA256}" "${ACQUISITION}/SHA256SUMS" \
  73d6a2083d763e9afd1eee2276a1df08119ebf9e40ebaf4115e95d32b59f44e8 "${CONFIG}/checkpoints.json" \
  142dcf3e57361e81e43f1fd8fd715d776a915ba43b4bd0dd047b76aa6b247eca "${CONFIG}/exposure_audit.json" \
  1688aa2b8a535bde4f3b3c2de11fc9fb0403117657e79913d61a99c719a24709 "${CONFIG}/development_crosswalk.json" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
CONTROL_PYTHON=$(command -v python)
PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - "${RUNTIME}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY

PYTHONPATH="${ROOT}:${ROOT}/src" "${CONTROL_PYTHON}" -m unittest -v \
  tests.contract.test_registry_files.RegistryFileContractTest.test_scimilarity_preflight_and_exposure_are_frozen \
  > "${STAGE}/validation/contract-test.stdout.txt" \
  2> "${STAGE}/validation/contract-test.stderr.txt"
tar -xzf "${RUNTIME}/payload/scimilarity-v1-1-compat-env.tar.gz" -C "${WORK}/runtime"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${WORK}/runtime/bin/python" -m unittest -v \
  tests.unit.test_scimilarity_v1_1_frozen_screen_50000 \
  > "${STAGE}/validation/unit-tests.stdout.txt" \
  2> "${STAGE}/validation/unit-tests.stderr.txt"

PYTHONPATH="${ROOT}:${ROOT}/src" "${CONTROL_PYTHON}" scripts/audit_scimilarity_v1_1_activation.py \
  --archive "${ACQUISITION}/model_v1.1.tar.gz" \
  --model-root "${MODEL}" \
  --code-license "${WORK}/runtime/lib/python3.10/site-packages/scimilarity-0.4.1.dist-info/licenses/LICENSE" \
  --nn-models "${WORK}/runtime/lib/python3.10/site-packages/scimilarity/nn_models.py" \
  --checkpoint-config "${CONFIG}/checkpoints.json" \
  --exposure-audit "${CONFIG}/exposure_audit.json" \
  --crosswalk "${CONFIG}/development_crosswalk.json" \
  --model-registry "${MODEL_REGISTRY}" \
  --fixture-root "${FIXTURE}/fixture" \
  --output "${STAGE}/activation_audit_v2.json" \
  > "${STAGE}/validation/activation-audit.stdout.json"

"${WORK}/runtime/bin/python" - "${FIXTURE}/fixture" \
  "${MODEL}/encoder.ckpt" "${STAGE}/cpu_preflight.json" \
  "${STAGE}/validation/five_fold_embeddings.npy" <<'PY'
from pathlib import Path
from hashlib import sha256
import json
import sys

import numpy as np

from scripts.extract_scimilarity_v1_1_frozen_screen_50000 import (
    _load_encoder,
    embed_batch,
    load_fixture,
)

fixture = load_fixture(Path(sys.argv[1]))
model, loading = _load_encoder(Path(sys.argv[2]), device="cpu")
indices = np.asarray(
    [int(np.flatnonzero(fixture["outer_folds"] == fold)[0]) for fold in range(5)],
    dtype=np.int64,
)
first = embed_batch(model, fixture["matrix"][indices], device="cpu")
second = embed_batch(model, fixture["matrix"][indices], device="cpu")
if not np.array_equal(first, second):
    raise SystemExit("SCimilarity five-fold CPU embedding repeat differs")
np.save(sys.argv[4], first, allow_pickle=False)
digest = sha256(Path(sys.argv[4]).read_bytes()).hexdigest()
receipt = {
    "schema_version": "masld-bench-scimilarity-v1-1-frozen-screen-cpu-preflight-v1",
    "status": "pass_cpu_adapter_preflight",
    "model_id": "scimilarity_v1_1",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "rows_in_frozen_fixture": len(fixture["rows"]),
    "representative_rows_embedded": len(indices),
    "representative_outer_folds": list(range(5)),
    "embedding_width": int(first.shape[1]),
    "five_fold_embeddings_sha256": digest,
    "deterministic_repeat_equal": True,
    "checkpoint_sha256": "a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f",
    "strict_weights_only_restore": True,
    "parameter_count": loading["parameter_count"],
    "aggregate_development_exposure": "encoder_seen",
    "study_gene_observability_shape": list(fixture["study_observability"].shape),
    "study_gene_observability_is_global": bool(
        np.all(fixture["study_observability"] == fixture["observed_mask"][None, :])
    ),
    "observed_genes": int(fixture["observed_mask"].sum()),
    "native_model_accepts_observability_mask": False,
    "released_reference_index_used": False,
    "raw_unrefined_embeddings_only": True,
    "full_50000_extraction_submitted": False,
    "full_50000_extraction_queued": False,
    "gpu_used": False,
    "downstream_head_fit": False,
    "common_head_fit": False,
    "evaluation_label_columns_read": [],
    "histology_read": False,
    "sealed_outcomes_read": False,
    "sealed_champion_eligible": False,
}
Path(sys.argv[3]).write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
PY

"${WORK}/runtime/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  scripts/audit_scimilarity_v1_1_activation.py \
  scripts/build_scimilarity_v1_1_frozen_screen_50000.py \
  scripts/extract_scimilarity_v1_1_frozen_screen_50000.py \
  tests/unit/test_scimilarity_v1_1_frozen_screen_50000.py \
  tests/contract/test_registry_files.py \
  slurm/preflight_scimilarity_v1_1_frozen_screen_50000_cpu.sbatch \
  "${CONFIG}/checkpoints.json" \
  "${CONFIG}/exposure_audit.json" \
  "${CONFIG}/development_crosswalk.json" \
  "${MODEL_REGISTRY}" \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="${ROOT}/src" "${CONTROL_PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((stage / "activation_audit_v2.json").read_text(encoding="utf-8"))
receipt = json.loads((stage / "cpu_preflight.json").read_text(encoding="utf-8"))
if (
    audit.get("schema_version") != "masld-bench-scimilarity-v1-1-activation-audit-v2"
    or audit.get("status") != "pass_exact_assets_terms_and_exposure"
    or audit.get("central_registry_reconciled") is not True
    or audit.get("aggregate_development_exposure") != "encoder_seen"
    or audit.get("study_gene_observability_shape") != [7, 28231]
    or audit.get("study_gene_observability_is_global") is not True
    or audit.get("evaluation_labels_read") is not False
    or audit.get("sealed_outcomes_read") is not False
    or audit.get("sealed_champion_eligible") is not False
    or receipt.get("status") != "pass_cpu_adapter_preflight"
    or receipt.get("rows_in_frozen_fixture") != 50_000
    or receipt.get("representative_outer_folds") != [0, 1, 2, 3, 4]
    or receipt.get("embedding_width") != 128
    or receipt.get("deterministic_repeat_equal") is not True
    or receipt.get("study_gene_observability_shape") != [7, 28231]
    or receipt.get("study_gene_observability_is_global") is not True
    or receipt.get("gpu_used") is not False
    or receipt.get("full_50000_extraction_submitted") is not False
    or receipt.get("full_50000_extraction_queued") is not False
    or receipt.get("downstream_head_fit") is not False
    or receipt.get("common_head_fit") is not False
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("SCimilarity 50k CPU preflight receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "scimilarity_v1_1_frozen_screen_cpu_preflight",
        "model_id": "scimilarity_v1_1",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "development_rows": 50_000,
        "aggregate_development_exposure": "encoder_seen",
        "study_gene_observability_shape": [7, 28231],
        "study_gene_observability_is_global": True,
        "gpu_used": False,
        "full_50000_extraction_submitted": False,
        "full_50000_extraction_queued": False,
        "downstream_head_fit": False,
        "common_head_fit": False,
        "evaluation_labels_read": False,
        "sealed_outcomes_read": False,
        "sealed_champion_eligible": False,
        "status": "passed",
    },
)
PY

trap - EXIT
rm -rf -- "${WORK}"
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
