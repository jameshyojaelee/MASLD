#!/usr/bin/env bash
#SBATCH --job-name=model-check-221
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
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
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_sequence_context_refit_dag.toml
ROW_UNIVERSE=${ROOT}/executions/model-check-220-21088696/contract/row_universe.tsv
TARGET=${ROOT}/executions/model-check-221-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE281364 refit-DAG attempt: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src

"${PYTHON}" -m unittest -v \
  tests.unit.test_prepare_gse281364_sequence_context_refit_dag \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/prepare_gse281364_sequence_context_refit_dag.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/contract" \
  > "${STAGE}/validation/dag.stdout.json" \
  2> "${STAGE}/validation/dag.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/prepare_gse281364_sequence_context_refit_dag.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --materialize-corgi-contexts \
  --output "${STAGE}/corgi_context_roster" \
  > "${STAGE}/validation/context.stdout.json" \
  2> "${STAGE}/validation/context.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/project_gse281364_long_range_embeddings.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --row-universe "${ROW_UNIVERSE}" \
  --model-roster caduceus,dnabert2,hyenadna \
  --output "${STAGE}/open_sequence_projection" \
  > "${STAGE}/validation/projection.stdout.json" \
  2> "${STAGE}/validation/projection.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/prepare_gse281364_sequence_context_refit_dag.py" \
  "${ROOT}/scripts/project_gse281364_long_range_embeddings.py" \
  "${ROOT}/tests/unit/test_prepare_gse281364_sequence_context_refit_dag.py" \
  "${ROOT}/slurm/prepare_gse281364_sequence_context_refit_dag_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


stage = Path(sys.argv[1]).resolve(strict=True)
dag = json.loads((stage / "contract/receipt.json").read_text(encoding="utf-8"))
contexts = json.loads(
    (stage / "corgi_context_roster/receipt.json").read_text(encoding="utf-8")
)
projection = json.loads(
    (stage / "open_sequence_projection/receipt.json").read_text(encoding="utf-8")
)
if (
    dag["status"]
    != "prepared_open_sequence_refit_context_restricted_no_trigger_pair"
    or dag["inventory_rows"] != 12
    or dag["stage_rows"] != 12
    or dag["current_trigger_compatible_pair_count"] != 0
    or dag["outcome_member_opened"]
    or dag["outcomes_read"]
    or dag["prediction_values_read"]
    or dag["model_fit"]
    or dag["stack_fit"]
    or dag["residual_correlation_calculated"]
    or dag["conditional_model_built_or_fit"]
    or contexts["contexts"] != 5
    or contexts["context_width"] != 2891
    or contexts["condition_matched"]
    or contexts["outcomes_read"]
    or contexts["model_forward_executed"]
    or projection["models"] != ["caduceus", "dnabert2", "hyenadna"]
    or projection["elements"] != 1033
    or projection["long_range_blocks"] != 239
    or projection["outcomes_read"]
    or projection["prediction_values_read"]
    or projection["downstream_head_fit"]
):
    raise SystemExit("GSE281364 sequence/context refit-DAG receipt differs")

freeze_tree(
    stage,
    {
        "artifact_class": "gse281364_sequence_context_refit_dag",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "inventory_count": 12,
        "stage_count": 12,
        "open_sequence_projection_models": ["caduceus", "dnabert2", "hyenadna"],
        "cross_cohort_context_count": 5,
        "context_condition_matched": False,
        "current_trigger_compatible_pair_count": 0,
        "outcomes_read": False,
        "prediction_values_read": False,
        "model_fit": False,
        "stack_fit": False,
        "conditional_model_built_or_fit": False,
        "status": "prepared_open_sequence_refit_context_restricted_no_trigger_pair",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 sequence/context refit DAG: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
