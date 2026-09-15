#!/usr/bin/env bash
#SBATCH --job-name=model-predict-083
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CONTRACT=${ROOT}/config/evaluation/gse260666_external_histology_prediction_v1.json
CONTRACT_SHA256=1d537d44169be2d6f4e0d4c5f3cbd3d1d7c57d5638e6261b864dd6f37506b009
TARGET=${ROOT}/executions/model-predict-083-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE260666 prediction-only artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CONTRACT_SHA256}" "${CONTRACT}" \
  | sha256sum --check --strict
bash -n slurm/predict_gse260666_external_histology_transfer_cpu.sbatch
"${PYTHON}" -m unittest -v \
  tests.unit.test_predict_gse260666_external_histology_transfer \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" scripts/predict_gse260666_external_histology_transfer.py \
  --benchmark-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --contract-sha256 "${CONTRACT_SHA256}" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/prediction_stdout.json" 2> "${STAGE}/prediction_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/predict_gse260666_external_histology_transfer.py \
  tests/unit/test_predict_gse260666_external_histology_transfer.py \
  slurm/predict_gse260666_external_histology_transfer_cpu.sbatch \
  "${CONTRACT}" \
  "${RUNTIME}/ARTIFACTS.json" \
  executions/model-fit-082-21092343/source_fit/ARTIFACTS.json \
  executions/model-activation-081-21092266/activation/ARTIFACTS.json \
  executions/model-data-068-21080339/fixture/model_input/ARTIFACTS.json \
  > "${STAGE}/source.sha256"
if rg -n 'evaluator_only|labels\.tsv' "${STAGE}/source.sha256" "${CONTRACT}"; then
  printf 'Prediction-only source unexpectedly names evaluator outcomes.\n' >&2
  exit 1
fi
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.contracts import PredictionBundle

root = Path(sys.argv[1]).resolve(strict=True)
predictions = root / "predictions"
receipt = json.loads((predictions / "prediction_receipt.json").read_text(encoding="utf-8"))
with (predictions / "prediction_bundle_index.tsv").open(encoding="utf-8", newline="") as handle:
    bundles = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt.get("status") != "passed_frozen_prediction_only_no_label_access"
    or receipt.get("participants") != 16
    or receipt.get("prediction_bundles") != 4
    or receipt.get("external_labels_read") is not False
    or receipt.get("query_fit_or_calibration_performed") is not False
    or receipt.get("external_metrics_calculated") is not False
    or len(bundles) != 4
):
    raise SystemExit("GSE260666 prediction-only verification differs")
for row in bundles:
    bundle_root = predictions / row["model_id"]
    verify_frozen_tree(bundle_root)
    bundle = PredictionBundle.load_json(bundle_root / "prediction_bundle.json")
    bundle.validate_artifacts(bundle_root)
    if (
        bundle.n_predictions != 16
        or bundle.metadata.get("external_labels_read") is not False
        or bundle.metadata.get("query_fit_or_calibration_performed") is not False
    ):
        raise SystemExit("GSE260666 PredictionBundle differs")
freeze_tree(
    predictions,
    {
        "artifact_class": "gse260666_external_histology_predictions",
        "participants": 16,
        "prediction_bundles": 4,
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
        "artifact_class": "gse260666_external_histology_prediction_campaign",
        "external_labels_read": False,
        "query_fit_or_calibration_performed": False,
        "metrics_calculated": False,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/predictions/ARTIFACTS.json"
