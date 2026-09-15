#!/usr/bin/env bash
#SBATCH --job-name=model-eval-934
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
OUTCOMES=${ROOT}/executions/model-data-061-21079623/activation
OUTCOMES_SHA=98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a
PREDICTIONS=${ROOT}/executions/model-predict-933-21133916/predictions
PREDICTIONS_SHA=9c615e35fd6ca00677d24cb48d226abd47167bb61a8a763a113a72efa4ba76a8
GATE=${ROOT}/config/evaluation/gse267145_nas_sum_transfer_promotion_gate.json
GATE_SHA=09b900c4ad8167a0d2cb4fa5ce49fdd6a06fc918bd22cbba313deb154bd42c51
TASK_SPEC=${ROOT}/config/evaluation/gse267145_nas_sum_transfer_task.toml
CANDIDATE=gene_median_ridge
OUT=${ROOT}/executions/model-eval-934-${SLURM_JOB_ID}
STAGE=${OUT}.staging
FAILED=${OUT}.failed
MODELS=(gene_median_pca_ridge gene_median_ridge per_array_rank_pca_ridge per_array_rank_ridge training_mean_nas_sum)

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${OUT}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite NAS evaluation artifacts.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/evaluator_only"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${OUTCOMES_SHA}" "${OUTCOMES}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${PREDICTIONS_SHA}" "${PREDICTIONS}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${GATE_SHA}" "${GATE}" | sha256sum --check --strict
bash -n slurm/evaluate_gse267145_nas_sum_transfer_cpu.sbatch
"${PYTHON}" - "${TASK_SPEC}" "${GATE_SHA}" <<'PY'
import re, sys
text = open(sys.argv[1], encoding="utf-8").read()
m = re.search(r'promotion_gate_config_sha256\s*=\s*"([0-9a-f]{64})"', text)
if m is None or m.group(1) != sys.argv[2]:
    raise SystemExit("TaskSpec does not pin this promotion gate")
print("promotion gate matches the TaskSpec pin")
PY
"${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse267145_nas_sum_transfer \
  > "${STAGE}/unit_tests.txt" 2>&1
PRED_ARGS=()
: > "${STAGE}/pre_label_prediction_verification.txt"
for m in "${MODELS[@]}"; do
  sha256sum "${PREDICTIONS}/${m}/ARTIFACTS.json" >> "${STAGE}/pre_label_prediction_verification.txt"
  PRED_ARGS+=(--prediction "${PREDICTIONS}/${m}")
  PRED_ARGS+=(--prediction-artifacts-sha256 "$(sha256sum "${PREDICTIONS}/${m}/ARTIFACTS.json" | cut -d' ' -f1)")
done
"${PYTHON}" scripts/build_gse267145_nas_sum_evaluator_labels.py \
  --endpoints "${OUTCOMES}/participant_endpoints.tsv" \
  --expected-participants 99 \
  $(for m in "${MODELS[@]}"; do printf -- '--prediction %s --prediction-artifacts-sha256 %s ' \
      "${PREDICTIONS}/${m}" "$(sha256sum "${PREDICTIONS}/${m}/ARTIFACTS.json" | cut -d' ' -f1)"; done) \
  --output-labels "${STAGE}/evaluator_only/labels.tsv" \
  --output-receipt "${STAGE}/evaluator_only/receipt.json" \
  > "${STAGE}/label_stdout.json" 2> "${STAGE}/label_stderr.txt"
"${PYTHON}" - "${STAGE}/evaluator_only" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "gse267145_nas_sum_evaluator_only_outcome",
    "participants": 99, "outcome_is_evaluator_only": True,
    "materialized_after_prediction_freeze": True, "status": "passed",
})
verify_frozen_tree(root)
PY
EVAL_SHA=$(sha256sum "${STAGE}/evaluator_only/ARTIFACTS.json" | cut -d' ' -f1)
"${PYTHON}" scripts/evaluate_gse267145_nas_sum_transfer.py \
  --evaluator-only "${STAGE}/evaluator_only" \
  --evaluator-only-artifacts-sha256 "${EVAL_SHA}" \
  --gate "${GATE}" --gate-sha256 "${GATE_SHA}" \
  --candidate-model-id "${CANDIDATE}" \
  "${PRED_ARGS[@]}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/evaluation_stdout.json" 2> "${STAGE}/evaluation_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/evaluate_gse267145_nas_sum_transfer.py \
  scripts/build_gse267145_nas_sum_evaluator_labels.py \
  tests/unit/test_evaluate_gse267145_nas_sum_transfer.py \
  slurm/evaluate_gse267145_nas_sum_transfer_cpu.sbatch \
  "${GATE}" "${TASK_SPEC}" "${RUNTIME}/ARTIFACTS.json" \
  "${OUTCOMES}/ARTIFACTS.json" "${PREDICTIONS}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
ev = root / "evaluation"
r = json.loads((ev / "receipt.json").read_text(encoding="utf-8"))
v = json.loads((ev / "promotion_gate_verdict.json").read_text(encoding="utf-8"))
if (
    r.get("status") != "passed_external_development_scoring"
    or r.get("participants") != 99
    or r.get("endpoint_binarised") is not False
    or r.get("recentred_or_rescaled_onto_the_target") is not False
    or r.get("model_selected_from_target_outcomes") is not False
    or v.get("verdict") not in {"PASS", "FAIL"}
    or v.get("conditions_total") != 4
    or v.get("conditions_applicable") > v.get("conditions_total")
):
    raise SystemExit("NAS evaluation verification differs")
freeze_tree(ev, {
    "artifact_class": "gse267145_nas_sum_transfer_evaluation", "participants": 99,
    "promotion_gate_verdict": v["verdict"],
    "promotion_gate_verdict_form": v["verdict_form"],
    "external_development_only": True, "status": "passed",
})
verify_frozen_tree(ev)
freeze_tree(root, {
    "artifact_class": "gse267145_nas_sum_transfer_evaluation_campaign",
    "promotion_gate_verdict": v["verdict"], "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps(v, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${OUT}"
sha256sum "${OUT}/ARTIFACTS.json" "${OUT}/evaluation/ARTIFACTS.json" "${OUT}/evaluator_only/ARTIFACTS.json"
column -t "${OUT}/evaluation/metrics.tsv"
