#!/usr/bin/env bash
#SBATCH --job-name=model-predict-933
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CONTRACT=${ROOT}/config/evaluation/gse267145_nas_sum_transfer_prediction_v1.json
CONTRACT_SHA=9b2ce760c8c4502cc8ccf41d99ebffd18778319a5e74c6eb0a76dffbb93c7e8d
SOURCE_FIT=${ROOT}/executions/model-fit-931-21133778/source_fit
ACTIVATION=${ROOT}/executions/model-data-929-21133561/activation
TARGET=${ROOT}/executions/model-data-064-21079902/fixture/molecular
OUT=${ROOT}/executions/model-predict-933-${SLURM_JOB_ID}
STAGE=${OUT}.staging
FAILED=${OUT}.failed
MODELS=(gene_median_pca_ridge gene_median_ridge per_array_rank_pca_ridge per_array_rank_ridge training_mean_nas_sum)

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${OUT}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite NAS prediction artifacts.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${CONTRACT_SHA}" "${CONTRACT}" | sha256sum --check --strict
bash -n slurm/predict_gse267145_nas_sum_transfer_cpu.sbatch
# The predictor may not READ the target histology table.  Naming the endpoint id
# is unavoidable and naming the outcome columns is how they are blocked, so the
# guard targets the paths an outcome could actually arrive through.
if grep -n 'participant_endpoints\|model-data-061\|evaluator_only\|labels\.tsv' \
    scripts/predict_gse267145_nas_sum_transfer.py "${CONTRACT}"; then
  printf 'Predictor unexpectedly reads from the target outcome table.\n' >&2; exit 1
fi
"${PYTHON}" -m unittest -v \
  tests.unit.test_predict_gse267145_nas_sum_transfer \
  tests.unit.test_fit_gse135251_nas_sum_transfer_source_models \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/predict_gse267145_nas_sum_transfer.py \
  --benchmark-root "${ROOT}" --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA}" --output "${STAGE}/predictions" \
  > "${STAGE}/prediction_stdout.json" 2> "${STAGE}/prediction_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/predict_gse267145_nas_sum_transfer.py \
  scripts/fit_gse135251_nas_sum_transfer_source_models.py \
  tests/unit/test_predict_gse267145_nas_sum_transfer.py \
  slurm/predict_gse267145_nas_sum_transfer_cpu.sbatch \
  "${CONTRACT}" "${RUNTIME}/ARTIFACTS.json" "${SOURCE_FIT}/ARTIFACTS.json" \
  "${ACTIVATION}/ARTIFACTS.json" "${TARGET}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv, json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.contracts import PredictionBundle
root = Path(sys.argv[1]).resolve(strict=True)
pred = root / "predictions"
r = json.loads((pred / "prediction_receipt.json").read_text(encoding="utf-8"))
with (pred / "prediction_bundle_index.tsv").open(encoding="utf-8", newline="") as h:
    bundles = list(csv.DictReader(h, delimiter="\t"))
if (
    r.get("status") != "passed_frozen_prediction_only_no_outcome_access"
    or r.get("participants") != 99 or r.get("shared_genes") != 38787
    or r.get("prediction_bundles") != 5 or len(bundles) != 5
    or r.get("target_outcomes_read") is not False
    or r.get("target_histology_table_read") is not False
    or r.get("recentred_or_rescaled_onto_the_target") is not False
    or r.get("endpoint_binarised") is not False
    or r.get("target_metrics_calculated") is not False
):
    raise SystemExit("NAS prediction verification differs")
for row in bundles:
    b = pred / row["model_id"]
    verify_frozen_tree(b)
    bundle = PredictionBundle.load_json(b / "prediction_bundle.json")
    bundle.validate_artifacts(b)
    if bundle.n_predictions != 99 or bundle.metadata.get("target_outcomes_read") is not False:
        raise SystemExit("NAS PredictionBundle differs")
    with (b / "predictions.tsv").open(encoding="utf-8", newline="") as h:
        table = list(csv.DictReader(h, delimiter="\t"))
    if {"nash_crn_component_sum", "steatosis", "ballooning"} & set(table[0]):
        raise SystemExit("prediction table carries the outcome")
freeze_tree(pred, {
    "artifact_class": "gse267145_nas_sum_predictions", "participants": 99,
    "prediction_bundles": 5, "target_outcomes_read": False,
    "metrics_calculated": False, "status": "passed_unscored",
})
verify_frozen_tree(pred)
freeze_tree(root, {
    "artifact_class": "gse267145_nas_sum_prediction_campaign",
    "target_outcomes_read": False, "metrics_calculated": False,
    "status": "passed_unscored",
})
verify_frozen_tree(root)
print(json.dumps({
    "min_pairwise_correlation_across_seeds": r["min_pairwise_correlation_across_seeds"],
    "predicted_range": r["predicted_range"],
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${OUT}"
sha256sum "${OUT}/ARTIFACTS.json" "${OUT}/predictions/ARTIFACTS.json"
for m in "${MODELS[@]}"; do sha256sum "${OUT}/predictions/${m}/ARTIFACTS.json"; done
