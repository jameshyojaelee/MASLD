#!/usr/bin/env bash
#SBATCH --job-name=model-scoring-078
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
: "${ROOT_VERIFIED_PRODUCTION_ARTIFACTS_SHA256:?root-verified production ARTIFACTS SHA-256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
PRODUCTION=${ROOT}/executions/model-training-069-production-v4
PRODUCTION_SHA256=ebbe02b9cf0e11ea6d01ccf981652631801b382fbd12225710b9fe17013da59d
PRODUCTION_AGGREGATE_SHA256=fa6b0ca7512c7ec1bd615d7f981b658d1d0c9bab2053ac6254b55cf8eb0458ac
SOURCE_VALIDATION=${ROOT}/executions/model-check-079-21092029
SOURCE_VALIDATION_SHA256=4c57942489e0bece01b8c6c233ab37b80f4911c01d76b18bad9ab4d3b32e5d95
BOUND_VALIDATION=${ROOT}/executions/model-check-079-21092115
BOUND_VALIDATION_SHA256=601aa3e128177a13c014f31bb07fe37dddad1a28bdc0a70f008b644393651291
OUTCOMES=${ROOT}/executions/model-data-061-21079623/activation
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
REFERENCE=${ROOT}/scripts/evaluate_gse267145_histology_predictions.py
TASK_CONTRACT=${ROOT}/config/evaluation/gse267145_histology_state_task.toml
SURFACE=${ROOT}/config/evaluation/gse267145_histology_benchmark_surface.json
SCORING_CONTRACT=${ROOT}/config/evaluation/gse267145_histology_production_v4_scoring.json
SCORING_CONTRACT_SHA256=532f66eab44561b4b84124c745b24a839f7d267c970fe5fb1fa9206513be7971
SCORER=${ROOT}/scripts/score_gse267145_histology_production_v4.py
SCORER_SHA256=dad1fa5b68686eaea6eeaa9c9725a4b0aed23dd4e96fcedbe1857d02ea070a43
AUDITOR=${ROOT}/scripts/audit_gse267145_histology_production_v4_scores.py
AUDITOR_SHA256=51c527b49fbd96a2a01c3d04a2353d977891b717cd0978f1e091e7203e340772
TARGET=${ROOT}/executions/model-scoring-078-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

if [[ "${SCORING_CONTRACT_SHA256}" == PENDING* \
   || "${SCORER_SHA256}" == PENDING* \
   || "${AUDITOR_SHA256}" == PENDING* ]]; then
  printf 'Scoring wrapper remains unbound and cannot be submitted.\n' >&2
  exit 1
fi
if [[ "${ROOT_VERIFIED_PRODUCTION_ARTIFACTS_SHA256}" != "${PRODUCTION_SHA256}" ]]; then
  printf 'Root-verified production SHA-256 differs from bound scoring wrapper.\n' >&2
  exit 1
fi
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 v4 production scores.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

for pair in \
  "${PRODUCTION_SHA256} ${PRODUCTION}/ARTIFACTS.json" \
  "${PRODUCTION_AGGREGATE_SHA256} ${PRODUCTION}/aggregate/ARTIFACTS.json" \
  "${SOURCE_VALIDATION_SHA256} ${SOURCE_VALIDATION}/ARTIFACTS.json" \
  "${BOUND_VALIDATION_SHA256} ${BOUND_VALIDATION}/ARTIFACTS.json" \
  "${SCORING_CONTRACT_SHA256} ${SCORING_CONTRACT}" \
  "${SCORER_SHA256} ${SCORER}" \
  "${AUDITOR_SHA256} ${AUDITOR}" \
  "1d753b8e23fda8f9582e4d8cc2ce9c90eeae65d370ce342bf9d32ffdd2ba5920 ${REFERENCE}" \
  "53a85a843e75b1c0fa088af86cb2f2593976d28039807dfa17d1ec571603778c ${TASK_CONTRACT}" \
  "23367a28bd16266c79adde854534b124cfad7d2062788625839e78f3267a8c0b ${SURFACE}" \
  "98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a ${OUTCOMES}/ARTIFACTS.json" \
  "9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3 ${FOLDS}/ARTIFACTS.json" \
  "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb ${RUNTIME}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

"${PYTHON}" - "${SOURCE_VALIDATION}" "${BOUND_VALIDATION}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY

"${PYTHON}" "${SCORER}" \
  --production "${PRODUCTION}" \
  --production-artifacts-sha256 "${PRODUCTION_SHA256}" \
  --outcomes "${OUTCOMES}" \
  --folds "${FOLDS}" \
  --reference-evaluator "${REFERENCE}" \
  --task-contract "${TASK_CONTRACT}" \
  --surface "${SURFACE}" \
  --scoring-contract "${SCORING_CONTRACT}" \
  --scoring-contract-sha256 "${SCORING_CONTRACT_SHA256}" \
  --output "${STAGE}/scores" \
  > "${STAGE}/validation/scorer_stdout.json" \
  2> "${STAGE}/validation/scorer_stderr.txt"
SCORES_SHA256=$(sha256sum "${STAGE}/scores/ARTIFACTS.json" | awk '{print $1}')

"${PYTHON}" "${AUDITOR}" \
  --production "${PRODUCTION}" \
  --production-artifacts-sha256 "${PRODUCTION_SHA256}" \
  --scores "${STAGE}/scores" \
  --scores-artifacts-sha256 "${SCORES_SHA256}" \
  --scorer-source "${SCORER}" \
  --scorer-source-sha256 "${SCORER_SHA256}" \
  --outcomes "${OUTCOMES}" \
  --folds "${FOLDS}" \
  --reference-evaluator "${REFERENCE}" \
  --task-contract "${TASK_CONTRACT}" \
  --surface "${SURFACE}" \
  --scoring-contract "${SCORING_CONTRACT}" \
  --scoring-contract-sha256 "${SCORING_CONTRACT_SHA256}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/auditor_stdout.json" \
  2> "${STAGE}/validation/auditor_stderr.txt"
AUDIT_SHA256=$(sha256sum "${STAGE}/audit/ARTIFACTS.json" | awk '{print $1}')

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  "${PRODUCTION}/ARTIFACTS.json" \
  "${OUTCOMES}/ARTIFACTS.json" \
  "${FOLDS}/ARTIFACTS.json" \
  "${REFERENCE}" \
  "${TASK_CONTRACT}" \
  "${SURFACE}" \
  "${SCORING_CONTRACT}" \
  "${SCORER}" \
  "${AUDITOR}" \
  slurm/run_gse267145_histology_production_v4_scores_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"

"${PYTHON}" - "${STAGE}" "${PRODUCTION_SHA256}" "${SCORES_SHA256}" "${AUDIT_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
production_sha256, scores_sha256, audit_sha256 = sys.argv[2:]
verify_frozen_tree(root / "scores")
verify_frozen_tree(root / "audit")
scores = json.loads((root / "scores/scoring_receipt.json").read_text(encoding="utf-8"))
audit = json.loads((root / "audit/audit_receipt.json").read_text(encoding="utf-8"))
if (
    scores.get("production_artifacts_sha256") != production_sha256
    or scores.get("standardized_prediction_bundles") != 55
    or scores.get("source_stage5_endpoint_dispositions") != 11
    or scores.get("standardized_metric_artifacts") != 55
    or scores.get("standardized_metric_rows") != 495
    or scores.get("endpoint_native_computed_metric_rows") != 99
    or scores.get("non_native_not_applicable_metric_rows") != 396
    or scores.get("seed_or_fold_results_emitted") is not False
    or scores.get("biological_resampling_unit") != "participant"
    or scores.get("champion_claim_allowed") is not False
    or scores.get("external_transfer_claim_allowed") is not False
    or scores.get("diagnostic_or_prognostic_claim_allowed") is not False
    or audit.get("scores_artifacts_sha256") != scores_sha256
    or audit.get("prediction_bundles_verified") != 55
    or audit.get("source_stage5_endpoint_dispositions_verified") != 11
    or audit.get("metric_rows_verified") != 495
    or audit.get("endpoint_native_computed_metric_rows_verified") != 99
    or audit.get("non_native_not_applicable_metric_rows_verified") != 396
    or audit.get("model_ranking_performed") is not False
    or audit.get("champion_claim_allowed") is not False
    or audit.get("external_transfer_claim_allowed") is not False
):
    raise SystemExit("frozen score or independent audit receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_production_v4_scoring_campaign",
        "production_artifacts_sha256": production_sha256,
        "scores_artifacts_sha256": scores_sha256,
        "audit_artifacts_sha256": audit_sha256,
        "prediction_bundles": 55,
        "source_stage5_endpoint_dispositions": 11,
        "metric_rows": 495,
        "endpoint_native_computed_metric_rows": 99,
        "non_native_not_applicable_metric_rows": 396,
        "biological_resampling_unit": "participant",
        "model_ranking_performed": False,
        "champion_claim_allowed": False,
        "external_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "status": "passed_development_scoring_and_independent_audit",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/scores/ARTIFACTS.json" \
  "${TARGET}/audit/ARTIFACTS.json"
