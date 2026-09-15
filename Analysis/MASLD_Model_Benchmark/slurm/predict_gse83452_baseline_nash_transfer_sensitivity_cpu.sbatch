#!/usr/bin/env bash
#SBATCH --job-name=model-predict-887
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
ARM_ID=nash_vs_nafl_plus_control
CONTRACT=${ROOT}/config/evaluation/gse83452_baseline_nash_transfer_prediction_nash_vs_nafl_plus_control_v1.json
CONTRACT_SHA256=6200d8b9eb9c0a00313e8ad1ee625ef2728d480002141dc3af5fee4d4be02b03
SOURCE_FIT=${ROOT}/executions/model-fit-883-21130286-nash_vs_nafl_plus_control/source_fit
ACTIVATION=${ROOT}/executions/model-data-881-21130281-activation/activation
SUMMARY=${ROOT}/executions/model-data-088-21112137-summarization/summary
TARGET=${ROOT}/executions/model-predict-887-${SLURM_JOB_ID}-nash_vs_nafl_plus_control
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE83452 prediction-only artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/predict_gse83452_baseline_nash_transfer_sensitivity_cpu.sbatch
printf '%s  %s\n' "${CONTRACT_SHA256}" "${CONTRACT}" | sha256sum --check --strict
# The predictor must never name the record table or the evaluator outcome.
if grep -n 'gse83452_records\.tsv\|evaluator_only\|labels\.tsv' \
    scripts/predict_gse83452_baseline_nash_transfer.py "${CONTRACT}"; then
  printf 'Prediction-only source unexpectedly names the outcome table.\n' >&2
  exit 1
fi
"${PYTHON}" -m unittest -v \
  tests.unit.test_predict_gse83452_baseline_nash_transfer \
  tests.unit.test_fit_gse135251_nash_transfer_source_models \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/predict_gse83452_baseline_nash_transfer.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA256}" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/prediction_stdout.json" 2> "${STAGE}/prediction_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/predict_gse83452_baseline_nash_transfer.py \
  scripts/fit_gse135251_nash_transfer_source_models.py \
  tests/unit/test_predict_gse83452_baseline_nash_transfer.py \
  slurm/predict_gse83452_baseline_nash_transfer_sensitivity_cpu.sbatch \
  "${CONTRACT}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE_FIT}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" \
  "${SUMMARY}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${ARM_ID}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.contracts import PredictionBundle

root = Path(sys.argv[1]).resolve(strict=True)
arm_id = sys.argv[2]
predictions = root / "predictions"
receipt = json.loads((predictions / "prediction_receipt.json").read_text(encoding="utf-8"))
with (predictions / "prediction_bundle_index.tsv").open(encoding="utf-8", newline="") as handle:
    bundles = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt.get("status") != "passed_frozen_prediction_only_no_label_access"
    or receipt.get("arm_id") != arm_id
    or receipt.get("participants") != 152
    or receipt.get("timepoint") != "baseline"
    or receipt.get("followup_arrays_read") is not False
    or receipt.get("shared_genes") != 21465
    or receipt.get("external_labels_read") is not False
    or receipt.get("external_participant_table_read") is not False
    or receipt.get("query_fit_or_calibration_performed") is not False
    or receipt.get("external_metrics_calculated") is not False
    or "age_sex_logistic" not in receipt.get("unfittable_model_ids", [])
    or len(bundles) != receipt["prediction_bundles"]
):
    raise SystemExit("GSE83452 prediction-only verification differs")
for row in bundles:
    bundle_root = predictions / row["model_id"]
    verify_frozen_tree(bundle_root)
    bundle = PredictionBundle.load_json(bundle_root / "prediction_bundle.json")
    bundle.validate_artifacts(bundle_root)
    if (
        bundle.n_predictions != 152
        or bundle.metadata.get("external_labels_read") is not False
        or bundle.metadata.get("query_fit_or_calibration_performed") is not False
        or bundle.metadata.get("prediction_frozen_before_evaluator_label_join") is not True
    ):
        raise SystemExit("GSE83452 PredictionBundle differs")
    with (bundle_root / "predictions.tsv").open(encoding="utf-8", newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    if {"nash_status", "age", "sex", "intervention"} & set(table[0]):
        raise SystemExit("prediction table carries an outcome or covariate column")
freeze_tree(
    predictions,
    {
        "artifact_class": "gse83452_baseline_nash_transfer_predictions",
        "arm_id": arm_id,
        "participants": 152,
        "prediction_bundles": len(bundles),
        "external_labels_read": False,
        "query_fit_or_calibration_performed": False,
        "metrics_calculated": False,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(predictions)
freeze_tree(
    root,
    {
        "artifact_class": "gse83452_baseline_nash_transfer_prediction_campaign",
        "arm_id": arm_id,
        "external_labels_read": False,
        "query_fit_or_calibration_performed": False,
        "metrics_calculated": False,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(root)
print(json.dumps({
    "models": receipt["models"],
    "distinct_seed_prediction_vectors": receipt["distinct_seed_prediction_vectors"],
    "representation_spread_diagnostics": receipt["representation_spread_diagnostics"],
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/predictions/ARTIFACTS.json"
for model in $(awk -F'\t' 'NR>1{print $1}' "${TARGET}/predictions/prediction_bundle_index.tsv"); do
  sha256sum "${TARGET}/predictions/${model}/ARTIFACTS.json"
done
