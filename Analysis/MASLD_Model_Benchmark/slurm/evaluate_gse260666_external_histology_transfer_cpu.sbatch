#!/usr/bin/env bash
#SBATCH --job-name=model-eval-084
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CONTRACT=${ROOT}/config/evaluation/gse260666_external_histology_evaluation_v1.json
CONTRACT_SHA256=7c6cc46aee1538054cb2c6d540d2ecf7a4610e7084431d555bb371d83c4b1189
TASK=${ROOT}/config/evaluation/gse260666_external_histology_transfer_task.toml
TASK_SHA256=dcd4a9e83c83d625d03c19b73de21185380aa66118fb51d108ce74069b813e40
PREDICTION=${ROOT}/executions/model-predict-083-21096324
PREDICTION_SHA256=e069c6b7d5975fb831ebe7f22d6e01e24a51b20718af2f9e4da111bcfdaf3df3
PREDICTIONS=${PREDICTION}/predictions
PREDICTIONS_SHA256=713352c97eb491160c39b5717619294d99333e7952b45ecf676e625725b9ce9b
EVALUATOR_ONLY=${ROOT}/executions/model-data-068-21080339/fixture/evaluator_only
EVALUATOR_ONLY_SHA256=6ed2b6166d00d16420454bc55a31ebd75b24d2dd3220d8c05482a8cd38cb6213
LABELS_SHA256=4cfc056da00a1ea9429222fe75a57f1287fe31f89656f6d19d1d6df426303fb6
TARGET=${ROOT}/executions/model-eval-084-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE260666 independent evaluation artifacts.\n' >&2
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
printf '%s  %s\n' "${TASK_SHA256}" "${TASK}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PREDICTION_SHA256}" "${PREDICTION}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PREDICTIONS_SHA256}" "${PREDICTIONS}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/evaluate_gse260666_external_histology_transfer_cpu.sbatch
"${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse260666_external_histology_transfer \
  > "${STAGE}/unit_tests.txt" 2>&1

# Verify every frozen prediction and its requirements before evaluator-only outcomes open.
"${PYTHON}" - "${ROOT}" "${CONTRACT}" "${STAGE}" <<'PY'
from pathlib import Path
import csv
import hashlib
import json
import sys

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.contracts import PredictionBundle

root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
stage = Path(sys.argv[3]).resolve(strict=True)
campaign_spec = contract["prediction_campaign"]
campaign = root / campaign_spec["path"]
predictions = campaign / "predictions"

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

if (
    contract.get("status") != "prediction_campaign_frozen_evaluator_ready"
    or campaign_spec.get("artifacts_sha256")
    != "e069c6b7d5975fb831ebe7f22d6e01e24a51b20718af2f9e4da111bcfdaf3df3"
    or campaign_spec.get("predictions_artifacts_sha256")
    != "713352c97eb491160c39b5717619294d99333e7952b45ecf676e625725b9ce9b"
    or campaign_spec.get("participants") != 16
    or campaign_spec.get("prediction_bundles") != 4
    or campaign_spec.get("external_labels_read") is not False
    or campaign_spec.get("query_fit_or_calibration_performed") is not False
    or campaign_spec.get("metrics_calculated") is not False
):
    raise SystemExit("frozen prediction campaign contract differs")
verify_frozen_tree(campaign)
verify_frozen_tree(predictions)
receipt = json.loads((predictions / "prediction_receipt.json").read_text(encoding="utf-8"))
with (predictions / "prediction_bundle_index.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    index = {row["model_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
bundles = contract["prediction_bundles"]
if len(bundles) != 4 or set(index) != {row["model_id"] for row in bundles}:
    raise SystemExit("prediction bundle roster differs")
for specification in bundles:
    bundle_root = root / specification["path"]
    verify_frozen_tree(bundle_root)
    bundle = PredictionBundle.load_json(bundle_root / "prediction_bundle.json")
    bundle.validate_artifacts(bundle_root)
    model_id = specification["model_id"]
    if (
        bundle.model_id != model_id
        or bundle.n_predictions != 16
        or specification.get("participant_rows") != 16
        or sha256_file(bundle_root / "ARTIFACTS.json")
        != specification["artifacts_sha256"]
        or sha256_file(bundle_root / "prediction_bundle.json")
        != specification["prediction_document_sha256"]
        or sha256_file(bundle_root / "predictions.tsv")
        != specification["prediction_table_sha256"]
        or index[model_id]["bundle_artifacts_sha256"]
        != specification["artifacts_sha256"]
        or index[model_id]["prediction_table_sha256"]
        != specification["prediction_table_sha256"]
    ):
        raise SystemExit(f"prediction bundle differs: {model_id}")
if (
    receipt.get("status") != "passed_frozen_prediction_only_no_label_access"
    or receipt.get("external_labels_read") is not False
    or receipt.get("query_fit_or_calibration_performed") is not False
    or receipt.get("external_metrics_calculated") is not False
):
    raise SystemExit("prediction-only receipt differs")
(stage / "pre_label_prediction_verification.json").write_text(
    json.dumps(
        {
            "status": "passed_all_predictions_frozen_before_evaluator_labels_opened",
            "prediction_campaign_artifacts_sha256": campaign_spec["artifacts_sha256"],
            "predictions_artifacts_sha256": campaign_spec["predictions_artifacts_sha256"],
            "prediction_bundles": 4,
            "participant_rows_per_bundle": 16,
            "total_prediction_rows": 64,
            "external_labels_read_during_verification": False,
        },
        indent=2,
        sort_keys=True,
    ) + "\n",
    encoding="utf-8",
)
PY

# Outcome access begins only after the pre-label prediction verification exists.
printf '%s  %s\n' "${EVALUATOR_ONLY_SHA256}" "${EVALUATOR_ONLY}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${LABELS_SHA256}" "${EVALUATOR_ONLY}/labels.tsv" \
  | sha256sum --check --strict
"${PYTHON}" - "${EVALUATOR_ONLY}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "passed_evaluator_only"
    or receipt.get("participants") != 16
    or receipt.get("mapping_is_one_to_one") is not True
    or receipt.get("nafl_nash_pooled") is not False
):
    raise SystemExit("evaluator-only authority differs")
PY

"${PYTHON}" scripts/evaluate_gse260666_external_histology_transfer.py \
  --evaluator-only "${EVALUATOR_ONLY}" \
  --evaluator-only-artifacts-sha256 "${EVALUATOR_ONLY_SHA256}" \
  --prediction "${PREDICTIONS}/training_stage_distribution" \
  --prediction-artifacts-sha256 bb9a8dc7e4ad37430205c4bd995bcfd56f6d9d3a24f8738853d5dc2385eec28f \
  --prediction "${PREDICTIONS}/rna_hvg_pca_nearest_centroid" \
  --prediction-artifacts-sha256 68d6748358f4f80c687c3c2f7444a95fdceedad53f0a6d8183d00a7ca8246844 \
  --prediction "${PREDICTIONS}/rna_hvg_pca_elastic_net" \
  --prediction-artifacts-sha256 081f6303bc78659a3a9967c768755e6b3287846cf4987bebce789f0b163d40e0 \
  --prediction "${PREDICTIONS}/rna_hvg_pca_linear_svm" \
  --prediction-artifacts-sha256 b2ecdd81bc29c437ad392f6874be404d8670664554c8d3151a4f0054304d7dc3 \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/evaluation_stdout.json" 2> "${STAGE}/evaluation_stderr.txt"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
evaluation = root / "evaluation"

def read_tsv(name: str) -> list[dict[str, str]]:
    with (evaluation / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

receipt = json.loads((evaluation / "receipt.json").read_text(encoding="utf-8"))
metrics = read_tsv("metrics.tsv")
intervals = read_tsv("bootstrap_intervals.tsv")
per_class = read_tsv("per_class.tsv")
class_intervals = read_tsv("per_class_f1_bootstrap_intervals.tsv")
models = {
    "training_stage_distribution",
    "rna_hvg_pca_nearest_centroid",
    "rna_hvg_pca_elastic_net",
    "rna_hvg_pca_linear_svm",
}
if (
    "Ran 2 tests" not in (root / "unit_tests.txt").read_text(encoding="utf-8")
    or (root / "evaluation_stderr.txt").read_text(encoding="utf-8")
    or receipt.get("status") != "passed_external_development_scoring"
    or receipt.get("participants") != 16
    or set(receipt.get("models", ())) != models
    or receipt.get("biological_unit") != "participant"
    or receipt.get("primary_metric") != "class_balanced_macro_f1"
    or receipt.get("bootstrap_method")
    != "source_class_stratified_participant_bootstrap"
    or receipt.get("bootstrap_replicates") != 10000
    or receipt.get("bootstrap_seed") != 20260824
    or receipt.get("nafl_nash_pooled") is not False
    or receipt.get("query_fit_or_calibration_performed") is not False
    or receipt.get("model_selected_from_external_outcomes") is not False
    or receipt.get("external_outcomes_used_for_model_repair") is not False
    or receipt.get("nominal_or_confirmatory_p_values_computed") is not False
    or receipt.get("project_sealed") is not False
    or receipt.get("external_development_only") is not True
    or receipt.get("champion_claim_eligible") is not False
    or receipt.get("clinical_claim_eligible") is not False
    or receipt.get("diagnostic_or_prognostic_claim_eligible") is not False
    or receipt.get("confirmatory_claim_eligible") is not False
    or len(metrics) != 24
    or len(intervals) != 24
    or len(per_class) != 12
    or len(class_intervals) != 12
    or {row["model_id"] for row in metrics} != models
    or (evaluation / "candidate_vs_strongest_baseline.tsv").exists()
):
    raise SystemExit("GSE260666 independent evaluation contract differs")
freeze_tree(
    evaluation,
    {
        "artifact_class": "gse260666_external_histology_independent_evaluation",
        "status": "passed_external_development_scoring",
        "participants": 16,
        "biological_unit": "participant",
        "models": 4,
        "metric_rows": 24,
        "participant_bootstrap_replicates": 10000,
        "project_sealed": False,
        "external_development_only": True,
        "champion_or_clinical_claim_eligible": False,
    },
)
verify_frozen_tree(evaluation)
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/evaluate_gse260666_external_histology_transfer.py \
  tests/unit/test_evaluate_gse260666_external_histology_transfer.py \
  slurm/evaluate_gse260666_external_histology_transfer_cpu.sbatch \
  "${CONTRACT}" \
  "${TASK}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${PREDICTION}/ARTIFACTS.json" \
  "${PREDICTIONS}/ARTIFACTS.json" \
  "${EVALUATOR_ONLY}/ARTIFACTS.json" \
  "${EVALUATOR_ONLY}/labels.tsv" \
  "${STAGE}/evaluation/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "gse260666_external_histology_independent_evaluation_campaign",
        "status": "passed_external_development_scoring",
        "prediction_campaign_verified_before_labels_opened": True,
        "participants": 16,
        "biological_unit": "participant",
        "models": 4,
        "metric_rows": 24,
        "participant_bootstrap_replicates": 10000,
        "project_sealed": False,
        "external_development_only": True,
        "champion_or_clinical_claim_eligible": False,
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/evaluation/ARTIFACTS.json"
