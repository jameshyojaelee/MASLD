#!/usr/bin/env bash
#SBATCH --job-name=model-eval-866
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
PARTICIPANTS=${ROOT}/executions/model-data-072-21080651-participants
PARTICIPANTS_ARTIFACTS_SHA256=6655f23b2092e163532113e74d75d3c2a93857ee870b4f0a94fdcb4fa7f9dbc0
PREDICTIONS=${ROOT}/executions/model-predict-865-21130175/predictions
PREDICTIONS_ARTIFACTS_SHA256=a7c5289edb3d5fd8c329e92315f20f9cc4d0e5b6a42ce5b28afaff74248bd81a
GATE=${ROOT}/config/evaluation/gse49541_fibrosis_transfer_promotion_gate.json
GATE_SHA256=791dfa3199c3bdbdeb74bfae429cc8820364862c7187f5d1dbe11b46c03dc384
TASK_SPEC=${ROOT}/config/evaluation/gse49541_fibrosis_transfer_task.toml
CANDIDATE=per_array_rank_elastic_net
SOURCE_FIT_RECEIPT=${ROOT}/executions/model-fit-869-21130053/source_fit/fit_receipt.json
ARM_ID=primary_f0f1_vs_f3f4_f2_excluded
TARGET=${ROOT}/executions/model-eval-866-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
MODELS=(gene_median_elastic_net gene_median_linear_svm per_array_rank_elastic_net training_prevalence)

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE49541 evaluation artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/evaluator_only"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PARTICIPANTS_ARTIFACTS_SHA256}" "${PARTICIPANTS}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${GATE_SHA256}" "${GATE}" | sha256sum --check --strict
printf '%s  %s\n' "${PREDICTIONS_ARTIFACTS_SHA256}" "${PREDICTIONS}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/evaluate_gse49541_fibrosis_transfer_cpu.sbatch
# The TaskSpec pins the check it is scored against; refuse a silent swap.
"${PYTHON}" - "${TASK_SPEC}" "${GATE_SHA256}" <<'PY'
import hashlib, re, sys
text = open(sys.argv[1], encoding="utf-8").read()
match = re.search(r'promotion_gate_config_sha256\s*=\s*"([0-9a-f]{64})"', text)
if match is None or match.group(1) != sys.argv[2]:
    raise SystemExit("TaskSpec does not pin this promotion gate")
print("promotion gate matches the TaskSpec pin")
PY
"${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse49541_fibrosis_transfer \
  > "${STAGE}/unit_tests.txt" 2>&1
# The only baseline that may be absent is the one the source fit recorded as
# unfittable, and the source fit is the authority on that, not this script.
"${PYTHON}" - "${SOURCE_FIT_RECEIPT}" <<'PY'
import json, sys
receipt = json.load(open(sys.argv[1], encoding="utf-8"))
if receipt["unfittable_model_ids"] != ["gene_median_pca_elastic_net"]:
    raise SystemExit(
        f"declared unfittable roster differs from the source fit: {receipt['unfittable_model_ids']}"
    )
if "per_array_rank_elastic_net" not in receipt["fitted_model_ids"]:
    raise SystemExit("the candidate model was not fitted on the source")
print(json.dumps({
    "unfittable_model_ids": receipt["unfittable_model_ids"],
    "selected_source_model_id": receipt["selected_source_model_id"],
    "no_source_model_exceeds_its_own_null": receipt["no_source_model_exceeds_its_own_null"],
}, sort_keys=True))
PY
# Record every prediction bundle's frozen digest BEFORE any outcome is read.
: > "${STAGE}/pre_label_prediction_verification.txt"
PREDICTION_ARGS=()
for model in "${MODELS[@]}"; do
  sha256sum "${PREDICTIONS}/${model}/ARTIFACTS.json" \
    >> "${STAGE}/pre_label_prediction_verification.txt"
  PREDICTION_ARGS+=(--prediction "${PREDICTIONS}/${model}")
  PREDICTION_ARGS+=(--prediction-artifacts-sha256 \
    "$(sha256sum "${PREDICTIONS}/${model}/ARTIFACTS.json" | cut -d' ' -f1)")
done
"${PYTHON}" scripts/build_gse49541_fibrosis_evaluator_labels.py \
  --participants "${PARTICIPANTS}/gse49541_participants.tsv" \
  $(for model in "${MODELS[@]}"; do printf -- '--prediction %s --prediction-artifacts-sha256 %s ' \
      "${PREDICTIONS}/${model}" \
      "$(sha256sum "${PREDICTIONS}/${model}/ARTIFACTS.json" | cut -d' ' -f1)"; done) \
  --output-labels "${STAGE}/evaluator_only/labels.tsv" \
  --output-receipt "${STAGE}/evaluator_only/receipt.json" \
  > "${STAGE}/label_stdout.json" 2> "${STAGE}/label_stderr.txt"
"${PYTHON}" - "${STAGE}/evaluator_only" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "gse49541_fibrosis_evaluator_only_outcome",
    "participants": 72,
    "outcome_is_evaluator_only": True,
    "materialized_after_prediction_freeze": True,
    "status": "passed",
})
verify_frozen_tree(root)
PY
EVALUATOR_ONLY_SHA256=$(sha256sum "${STAGE}/evaluator_only/ARTIFACTS.json" | cut -d' ' -f1)
"${PYTHON}" scripts/evaluate_gse49541_fibrosis_transfer.py \
  --evaluator-only "${STAGE}/evaluator_only" \
  --evaluator-only-artifacts-sha256 "${EVALUATOR_ONLY_SHA256}" \
  --gate "${GATE}" \
  --gate-sha256 "${GATE_SHA256}" \
  --candidate-model-id "${CANDIDATE}" \
  --arm-id "${ARM_ID}" \
  --gate-eligible-arm true \
  --unfittable-baseline gene_median_pca_elastic_net \
  "${PREDICTION_ARGS[@]}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/evaluation_stdout.json" 2> "${STAGE}/evaluation_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/evaluate_gse49541_fibrosis_transfer.py \
  scripts/build_gse49541_fibrosis_evaluator_labels.py \
  tests/unit/test_evaluate_gse49541_fibrosis_transfer.py \
  slurm/evaluate_gse49541_fibrosis_transfer_cpu.sbatch \
  "${GATE}" \
  "${TASK_SPEC}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${PARTICIPANTS}/ARTIFACTS.json" \
  "${PREDICTIONS}/ARTIFACTS.json" \
  "${SOURCE_FIT_RECEIPT}" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
evaluation = root / "evaluation"
receipt = json.loads((evaluation / "receipt.json").read_text(encoding="utf-8"))
verdict = json.loads((evaluation / "promotion_gate_verdict.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "passed_external_development_scoring"
    or receipt.get("participants") != 72
    or receipt.get("external_development_only") is not True
    or receipt.get("champion_claim_eligible") is not False
    or receipt.get("model_selected_from_external_outcomes") is not False
    or verdict.get("verdict") not in {"PASS", "FAIL"}
    or verdict.get("conditions_total") != 4
):
    raise SystemExit("GSE49541 evaluation verification differs")
freeze_tree(evaluation, {
    "artifact_class": "gse49541_fibrosis_transfer_evaluation",
    "participants": 72,
    "promotion_gate_verdict": verdict["verdict"],
    "external_development_only": True,
    "champion_claim_eligible": False,
    "status": "passed",
})
verify_frozen_tree(evaluation)
freeze_tree(root, {
    "artifact_class": "gse49541_fibrosis_transfer_evaluation_campaign",
    "promotion_gate_verdict": verdict["verdict"],
    "external_development_only": True,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps(verdict, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/evaluation/ARTIFACTS.json" "${TARGET}/evaluator_only/ARTIFACTS.json"
column -t "${TARGET}/evaluation/metrics.tsv"
