#!/usr/bin/env bash
#SBATCH --job-name=model-check-079
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
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/model-check-079-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 v4 score validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

for pair in \
  "dad1fa5b68686eaea6eeaa9c9725a4b0aed23dd4e96fcedbe1857d02ea070a43 scripts/score_gse267145_histology_production_v4.py" \
  "51c527b49fbd96a2a01c3d04a2353d977891b717cd0978f1e091e7203e340772 scripts/audit_gse267145_histology_production_v4_scores.py" \
  "77a116614ff124054c192dd808131c270509073eceb453fcdb2bb7b9b8797a41 tests/unit/test_score_gse267145_histology_production_v4.py" \
  "187e9ed663843f153c553e48abd1072fcc5e462c7936f7819f045e16a95f4f12 tests/unit/test_audit_gse267145_histology_production_v4_scores.py" \
  "532f66eab44561b4b84124c745b24a839f7d267c970fe5fb1fa9206513be7971 config/evaluation/gse267145_histology_production_v4_scoring.json" \
  "1d753b8e23fda8f9582e4d8cc2ce9c90eeae65d370ce342bf9d32ffdd2ba5920 scripts/evaluate_gse267145_histology_predictions.py" \
  "53a85a843e75b1c0fa088af86cb2f2593976d28039807dfa17d1ec571603778c config/evaluation/gse267145_histology_state_task.toml" \
  "23367a28bd16266c79adde854534b124cfad7d2062788625839e78f3267a8c0b config/evaluation/gse267145_histology_benchmark_surface.json" \
  "98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a executions/model-data-061-21079623/activation/ARTIFACTS.json" \
  "9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3 executions/model-data-064-21079902/fixture/folds/ARTIFACTS.json" \
  "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json" \
  "ebbe02b9cf0e11ea6d01ccf981652631801b382fbd12225710b9fe17013da59d executions/model-training-069-production-v4/ARTIFACTS.json" \
  "fa6b0ca7512c7ec1bd615d7f981b658d1d0c9bab2053ac6254b55cf8eb0458ac executions/model-training-069-production-v4/aggregate/ARTIFACTS.json" \
  "4c57942489e0bece01b8c6c233ab37b80f4911c01d76b18bad9ab4d3b32e5d95 executions/model-check-079-21092029/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

"${PYTHON}" -m py_compile \
  scripts/score_gse267145_histology_production_v4.py \
  scripts/audit_gse267145_histology_production_v4_scores.py
"${PYTHON}" -m unittest -v \
  tests.unit.test_score_gse267145_histology_production_v4 \
  tests.unit.test_audit_gse267145_histology_production_v4_scores \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/score_gse267145_histology_production_v4.py \
  scripts/audit_gse267145_histology_production_v4_scores.py \
  tests/unit/test_score_gse267145_histology_production_v4.py \
  tests/unit/test_audit_gse267145_histology_production_v4_scores.py \
  config/evaluation/gse267145_histology_production_v4_scoring.json \
  scripts/evaluate_gse267145_histology_predictions.py \
  config/evaluation/gse267145_histology_state_task.toml \
  config/evaluation/gse267145_histology_benchmark_surface.json \
  slurm/validate_gse267145_histology_production_v4_scores_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
import scripts.audit_gse267145_histology_production_v4_scores as auditor
import scripts.score_gse267145_histology_production_v4 as scorer

root = Path(sys.argv[1]).resolve(strict=True)
contract = json.loads(
    Path("config/evaluation/gse267145_histology_production_v4_scoring.json").read_text(
        encoding="utf-8"
    )
)
if (
    contract.get("status") != "bound_root_verified_production_artifacts"
    or contract.get("production_campaign", {}).get("artifacts_sha256")
    != "ebbe02b9cf0e11ea6d01ccf981652631801b382fbd12225710b9fe17013da59d"
    or contract.get("production_campaign", {}).get("aggregate_artifacts_sha256")
    != "fa6b0ca7512c7ec1bd615d7f981b658d1d0c9bab2053ac6254b55cf8eb0458ac"
    or contract.get("independent_source_validation", {}).get("artifacts_sha256")
    != "4c57942489e0bece01b8c6c233ab37b80f4911c01d76b18bad9ab4d3b32e5d95"
    or len(scorer._endpoint_bundles()) != 55
    or len(auditor._prediction_sets()) != 55
    or len(scorer.PREDICTED_ENDPOINTS) != 5
    or sum(len(value) for value in scorer.ENDPOINT_METRICS.values()) != 9
    or auditor.EXPECTED_PREDICTION_BUNDLES != 55
    or auditor.EXPECTED_STAGE5_DISPOSITIONS != 11
    or auditor.EXPECTED_METRIC_ROWS != 495
    or auditor.EXPECTED_COMPUTED_METRIC_ROWS != 99
    or auditor.EXPECTED_NOT_APPLICABLE_METRIC_ROWS != 396
):
    raise SystemExit("v4 scoring validation mapping or submission gate differs")
verify_frozen_tree(Path("executions/model-training-069-production-v4"))
verify_frozen_tree(Path("executions/model-check-079-21092029"))
receipt = {
    "schema_version": "masld-bench-gse267145-production-v4-score-validation-v1",
    "status": "passed_unscored_source_and_unit_validation",
    "scoring_contract_status": contract["status"],
    "production_artifacts_bound": True,
    "production_artifacts_sha256": (
        "ebbe02b9cf0e11ea6d01ccf981652631801b382fbd12225710b9fe17013da59d"
    ),
    "production_aggregate_artifacts_sha256": (
        "fa6b0ca7512c7ec1bd615d7f981b658d1d0c9bab2053ac6254b55cf8eb0458ac"
    ),
    "source_validation_rebind_required": False,
    "independent_source_validation_artifacts_sha256": (
        "4c57942489e0bece01b8c6c233ab37b80f4911c01d76b18bad9ab4d3b32e5d95"
    ),
    "production_scores_calculated": False,
    "unit_tests_passed": True,
    "syntax_validation_passed": True,
    "prediction_bundles_expected": 55,
    "source_stage5_endpoint_dispositions_expected": 11,
    "metric_rows_expected": 495,
    "endpoint_native_computed_metric_rows_expected": 99,
    "non_native_not_applicable_metric_rows_expected": 396,
    "seed_or_fold_results_emitted": False,
    "biological_resampling_unit": "participant",
    "model_ranking_allowed": False,
    "champion_claim_allowed": False,
    "external_transfer_claim_allowed": False,
    "diagnostic_or_prognostic_claim_allowed": False,
}
(root / "validation/receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_production_v4_score_validation",
        "production_artifacts_bound": True,
        "production_artifacts_sha256": (
            "ebbe02b9cf0e11ea6d01ccf981652631801b382fbd12225710b9fe17013da59d"
        ),
        "production_aggregate_artifacts_sha256": (
            "fa6b0ca7512c7ec1bd615d7f981b658d1d0c9bab2053ac6254b55cf8eb0458ac"
        ),
        "source_validation_rebind_required": False,
        "independent_source_validation_artifacts_sha256": (
            "4c57942489e0bece01b8c6c233ab37b80f4911c01d76b18bad9ab4d3b32e5d95"
        ),
        "production_scores_calculated": False,
        "prediction_bundles_expected": 55,
        "source_stage5_endpoint_dispositions_expected": 11,
        "metric_rows_expected": 495,
        "endpoint_native_computed_metric_rows_expected": 99,
        "non_native_not_applicable_metric_rows_expected": 396,
        "status": "passed_unscored",
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
