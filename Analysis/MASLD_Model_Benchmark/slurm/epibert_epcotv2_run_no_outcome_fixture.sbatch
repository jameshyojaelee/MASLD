#!/usr/bin/env bash
#SBATCH --job-name=masld-observed-atac-fixture
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/epibert-tensorflow-2.13.0-21064371"
RUNTIME_ARTIFACTS_SHA256=10a8691e9468dfa44a25e55770b477e77d2ebed12839d49203f1a386596af27f
SOURCES="${ROOT}/executions/epibert-epcotv2-admission-sources-21064417"
SOURCES_ARTIFACTS_SHA256=26e12915b05c2727c45eca091f85ec89f96c205e817c6ec9b64f181a563e7577
IMAGE="${RUNTIME}/tensorflow-tensorflow-2.13.0-gpu.sif"
TARGET="${ROOT}/executions/epibert-epcotv2-no-outcome-fixture-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite observed-ATAC fixture: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/fixtures" "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCES_ARTIFACTS_SHA256}" "${SOURCES}/ARTIFACTS.json" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${SOURCES}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
module purge
module load singularity/4.2.2
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/scripts/epibert_epcotv2_no_outcome_fixture.py" \
    --model epibert \
    --runtime tensorflow \
    --output "${STAGE}/fixtures/epibert" \
    > "${STAGE}/validation/epibert.stdout.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/scripts/epibert_epcotv2_no_outcome_fixture.py" \
    --model epcotv2 \
    --runtime numpy \
    --output "${STAGE}/fixtures/epcotv2" \
    > "${STAGE}/validation/epcotv2.stdout.txt"
singularity exec --cleanenv "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${ROOT}/scripts/epibert_epcotv2_no_outcome_fixture.py" \
  "${ROOT}/slurm/epibert_epcotv2_run_no_outcome_fixture.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" "${SOURCES_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
epibert = json.loads((root / "fixtures/epibert/receipt.json").read_text())
epcot = json.loads((root / "fixtures/epcotv2/receipt.json").read_text())
if (
    epibert["status"] != "pass"
    or epibert["runtime_input_tensor_probe"]["version"] != "2.13.0"
    or not epibert["runtime_input_tensor_probe"]["exact_native_version"]
    or epibert["runtime_input_tensor_probe"]["device"] != "/GPU:0"
    or epibert["outcomes_read"]
    or epibert["model_checkpoint_loaded"]
    or epibert["query_semantics"]["rna_conditioning"]
    or not epibert["query_semantics"]["scored_atac_span_removed_from_input"]
    or epcot["status"] != "pass"
    or epcot["outcomes_read"]
    or epcot["model_checkpoint_loaded"]
    or epcot["query_semantics"]["rna_conditioning"]
    or epcot["query_semantics"]["atac_is_output_head"]
):
    raise SystemExit("observed-ATAC no-outcome fixture differs")
freeze_tree(
    root,
    {
        "artifact_class": "observed_atac_no_outcome_runtime_fixture",
        "models": ["epibert", "epcotv2"],
        "epibert_terminal_disposition": "input_runtime_passed_checkpoint_restore_and_numeric_parity_pending",
        "epcotv2_terminal_disposition": "input_schema_passed_execution_blocked_terms",
        "runtime_artifacts_sha256": sys.argv[2],
        "admission_sources_artifacts_sha256": sys.argv[3],
        "outcomes_read": False,
        "model_checkpoints_loaded": False,
        "rna_conditioned_atac_models": [],
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Observed-ATAC no-outcome fixture: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
