#!/usr/bin/env bash
#SBATCH --job-name=model-preflight-507
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
CONTRACT=${ROOT}/config/artifacts/models/uce/frozen_screen_50000.json
CONTRACT_SHA256=f771b97b17a8b0be7f8657f1dd3943f0e605844d026ff1101c459edcaa0b4087
COUNT_ROOT=${ROOT}/executions/model-data-051-21077382
COUNT_ROOT_SHA256=ba40ce9555bb6d7bc278be2a80015b993f9f3650fa06f59505ad279c88cc117a
SYMBOL_ROOT=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SYMBOL_ROOT_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT_ROOT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_ROOT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-preflight-507-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${exit_code}"
}
trap finish EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite UCE preflight output.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${CONTRACT_SHA256}" "${CONTRACT}" \
  "${COUNT_ROOT_SHA256}" "${COUNT_ROOT}/ARTIFACTS.json" \
  "${SYMBOL_ROOT_SHA256}" "${SYMBOL_ROOT}/ARTIFACTS.json" \
  "${SPLIT_ROOT_SHA256}" "${SPLIT_ROOT}/ARTIFACTS.json" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"

export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" - "${COUNT_ROOT}" "${SYMBOL_ROOT}" "${SPLIT_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for path in map(Path, sys.argv[1:]):
    verify_frozen_tree(path)
PY
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_uce_frozen_screen_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${RUNTIME}/env/bin/python" scripts/build_uce_frozen_screen_50000.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/fixture" \
  --inspect-checkpoints \
  > "${STAGE}/validation/builder.stdout.txt" \
  2> "${STAGE}/validation/builder.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  config/artifacts/models/uce/checkpoints.json \
  config/artifacts/models/uce/exposure_audit.json \
  config/artifacts/models/uce/frozen_screen_50000.json \
  scripts/uce_frozen_screen_50000_common.py \
  scripts/build_uce_frozen_screen_50000.py \
  scripts/extract_uce_frozen_screen_50000.py \
  tests/unit/test_uce_frozen_screen_50000.py \
  slurm/preflight_uce_frozen_screen_50000_cpu.sbatch \
  slurm/run_uce_frozen_screen_50000.sbatch \
  > "${STAGE}/validation/source.sha256"

"${RUNTIME}/env/bin/python" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "fixture/receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_outcome_blind_uce_fixture"
    or receipt.get("rows") != 50_000
    or receipt.get("donors") != 102
    or receipt.get("studies") != 7
    or receipt.get("cells_removed") != 0
    or receipt.get("cells_truncated") != 0
    or set(receipt.get("checkpoint_inspection", {})) != {"uce_4l", "uce_33l"}
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("symbol_source_observation_fields_read")
    or receipt.get("sealed_outcomes_read") is not False
    or receipt.get("downstream_head_fit") is not False
):
    raise SystemExit("UCE fixture receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "uce_frozen_screen_50000_fixture",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "rows": 50_000,
        "donors": 102,
        "studies": 7,
        "evaluation_labels_read": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "downstream_head_fit": False,
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
