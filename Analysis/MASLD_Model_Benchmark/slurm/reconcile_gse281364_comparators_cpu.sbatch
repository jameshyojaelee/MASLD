#!/usr/bin/env bash
#SBATCH --job-name=model-check-225
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
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
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_comparator_reconciliation.toml
TARGET=${ROOT}/executions/model-check-225-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE281364 comparator reconciliation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src

"${PYTHON}" -m py_compile \
  scripts/reconcile_gse281364_comparators.py \
  tests/unit/test_reconcile_gse281364_comparators.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_reconcile_gse281364_comparators \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" scripts/reconcile_gse281364_comparators.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/reconciliation" \
  > "${STAGE}/validation/reconciliation.stdout.json" \
  2> "${STAGE}/validation/reconciliation.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/gse281364_comparator_reconciliation.toml \
  scripts/reconcile_gse281364_comparators.py \
  tests/unit/test_reconcile_gse281364_comparators.py \
  slurm/reconcile_gse281364_comparators_cpu.sbatch \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import csv
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


stage = Path(sys.argv[1]).resolve(strict=True)
root = stage / "reconciliation"
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
interpretation = json.loads(
    (root / "terminal_interpretation_audit.json").read_text(encoding="utf-8")
)
with (root / "comparator_disposition.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    dispositions = list(csv.DictReader(handle, delimiter="\t"))
with (root / "mpralegnet_selected_elements.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    selected_rows = sum(1 for _ in handle) - 1
with (root / "mpralegnet_row_key_join.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    joined_rows = sum(1 for _ in handle) - 1
if (
    receipt["status"] != "pass_outcome_blind_join_coverage_disposition"
    or selected_rows != 1033
    or joined_rows != 10330
    or receipt["long_range_blocks"] != 239
    or receipt["mpralegnet_reassigned_from_original_fold"] != 830
    or len(dispositions) != 5
    or dispositions[0]["model_id"] != "mpralegnet"
    or dispositions[0]["joined_row_keys"] != "10330"
    or any(row["valid_for_current_shortlist"] != "false" for row in dispositions)
    or receipt["candidate_summary_metric_values_read"]
    or receipt["metrics_calculated"]
    or receipt["outcomes_read"]
    or receipt["model_fit"]
    or receipt["new_predictions_generated"]
    or receipt["sealed_assets_read"]
    or not receipt["shortlist_blocked"]
    or not receipt["finalist_claim_blocked"]
    or not receipt["complementarity_blocked"]
    or not receipt["conditional_trigger_blocked"]
    or receipt["champion_claim"]
    or interpretation["candidate_summary_metric_values_read"]
    or interpretation[
        "positive_gain_seeds_for_deterministic_heads_is_stability_evidence"
    ]
    or interpretation[
        "four_of_five_stability_gate_satisfied_by_deterministic_repeats"
    ]
):
    raise SystemExit("GSE281364 comparator reconciliation receipt differs")
freeze_tree(stage, {
    "artifact_class": "gse281364_outcome_blind_comparator_reconciliation",
    "dataset_id": "gse281364",
    "task_id": "variant_to_regulation",
    "selected_elements": 1033,
    "row_keys": 10330,
    "long_range_blocks": 239,
    "mpralegnet_join_coverage": 1.0,
    "mpralegnet_rekey_complete": True,
    "sequence_cnn_gse281364_predictions": 0,
    "sequence_transformer_gse281364_predictions": 0,
    "gkmsvm_gse281364_predictions": 0,
    "deltasvm_gse281364_predictions": 0,
    "candidate_summary_metric_values_read": False,
    "deterministic_seed_repeats_are_stability_evidence": False,
    "metrics_calculated": False,
    "outcomes_read": False,
    "model_fit": False,
    "new_predictions_generated": False,
    "shortlist_blocked": True,
    "finalist_claim_blocked": True,
    "complementarity_blocked": True,
    "conditional_trigger_blocked": True,
    "status": "pass_outcome_blind_join_coverage_disposition",
})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 comparator reconciliation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
