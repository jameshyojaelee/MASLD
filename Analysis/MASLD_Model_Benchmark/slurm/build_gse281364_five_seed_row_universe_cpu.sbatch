#!/usr/bin/env bash
#SBATCH --job-name=model-check-220
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_five_seed_row_universe.toml
TARGET=${ROOT}/executions/model-check-220-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE281364 row-universe attempt: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src

"${PYTHON}" -m unittest -v \
  tests.unit.test_build_gse281364_five_seed_row_universe \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/build_gse281364_five_seed_row_universe.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/contract" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/build_gse281364_five_seed_row_universe.py" \
  "${ROOT}/tests/unit/test_build_gse281364_five_seed_row_universe.py" \
  "${ROOT}/slurm/build_gse281364_five_seed_row_universe_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "contract/receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "row_universe_defined_stacking_blocked"
    or receipt["source_locus_groups"] != 1033
    or receipt["selected_elements"] != 1033
    or receipt["base_rows_per_seed"] != 2066
    or receipt["fixed_seeds"] != [1103, 2909, 4721, 6673, 8111]
    or receipt["seeded_rows"] != 10330
    or receipt["long_range_blocks"] != 239
    or receipt["source_groups_reassigned_from_original_6kb_fold_map"] != 830
    or receipt["current_open_sequence_source_rows"] != 0
    or receipt["current_open_context_source_rows"] != 0
    or receipt["current_trigger_compatible_pair_rows"] != 0
    or receipt["outcomes_read"]
    or receipt["prediction_values_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_assets_read"]
    or receipt["predictions_generated"]
    or receipt["stack_fit"]
    or receipt["residual_correlation_calculated"]
    or receipt["conditional_model_built_or_fit"]
):
    raise SystemExit("GSE281364 row-universe receipt differs")

freeze_tree(
    stage,
    {
        "artifact_class": "gse281364_five_seed_row_universe",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "source_locus_group_count": 1033,
        "base_row_count_per_seed": 2066,
        "seeded_row_count": 10330,
        "long_range_block_count": 239,
        "fixed_seed_count": 5,
        "current_trigger_compatible_pair_count": 0,
        "outcomes_read": False,
        "prediction_values_read": False,
        "sealed_assets_read": False,
        "stack_fit": False,
        "conditional_model_built_or_fit": False,
        "status": "row_universe_defined_stacking_blocked",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 five-seed row universe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
