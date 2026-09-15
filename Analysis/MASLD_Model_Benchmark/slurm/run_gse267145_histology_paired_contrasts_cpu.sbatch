#!/usr/bin/env bash
#SBATCH --job-name=model-check-1102
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
: "${SEALED_PRESPEC:?sealed prespecification execution path is required}"
: "${SEALED_PRESPEC_ARTIFACTS_SHA256:?sealed prespecification ARTIFACTS SHA-256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
PRESPEC=${ROOT}/config/evaluation/gse267145_histology_paired_contrast_prespec_v1.json
PRESPEC_SHA256=fd00b3d2315fa9f3e729fb81e58ee3866c1f581c67d1aaed6a79c3b5cc7e094a
ANALYZER=${ROOT}/scripts/analyze_gse267145_histology_paired_contrasts.py
SCORES=${ROOT}/executions/model-scoring-078-21092130/scores
OUTCOMES=${ROOT}/executions/model-data-061-21079623/activation
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds
TARGET=${ROOT}/executions/model-check-1102-${SLURM_JOB_ID}-paired-contrasts
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the GSE267145 paired development contrasts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

for pair in \
  "${PRESPEC_SHA256} ${PRESPEC}" \
  "${SEALED_PRESPEC_ARTIFACTS_SHA256} ${SEALED_PRESPEC}/ARTIFACTS.json" \
  "dad1fa5b68686eaea6eeaa9c9725a4b0aed23dd4e96fcedbe1857d02ea070a43 ${ROOT}/scripts/score_gse267145_histology_production_v4.py" \
  "1d753b8e23fda8f9582e4d8cc2ce9c90eeae65d370ce342bf9d32ffdd2ba5920 ${ROOT}/scripts/evaluate_gse267145_histology_predictions.py" \
  "a55aed3d9ac4cad2211ed3b319c45f911d05851b8559408ed99bf7fa59562393 ${SCORES}/ARTIFACTS.json" \
  "98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a ${OUTCOMES}/ARTIFACTS.json" \
  "9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3 ${FOLDS}/ARTIFACTS.json" \
  "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb ${RUNTIME}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

"${PYTHON}" - "${SEALED_PRESPEC}" "${PRESPEC_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

sealed = Path(sys.argv[1]).resolve(strict=True)
prespec_sha256 = sys.argv[2]
verify_frozen_tree(sealed)
verify_frozen_tree(sealed / "prespec")
receipt = json.loads(
    (sealed / "prespec/prespec_seal_receipt.json").read_text(encoding="utf-8")
)
if (
    receipt.get("prespec_contract_sha256") != prespec_sha256
    or receipt.get("outcomes_read") is not False
    or receipt.get("contrasts_calculated") is not False
    or receipt.get("status") != "sealed_prespec_no_outcome_read"
):
    raise SystemExit("sealed prespecification does not bind this analysis")
PY

"${PYTHON}" "${ANALYZER}" \
  --prespec "${PRESPEC}" \
  --prespec-sha256 "${PRESPEC_SHA256}" \
  --sealed-prespec "${SEALED_PRESPEC}" \
  --sealed-prespec-artifacts-sha256 "${SEALED_PRESPEC_ARTIFACTS_SHA256}" \
  --scores "${SCORES}" \
  --outcomes "${OUTCOMES}" \
  --folds "${FOLDS}" \
  --output "${STAGE}/contrasts" \
  > "${STAGE}/validation/analyzer_stdout.json" \
  2> "${STAGE}/validation/analyzer_stderr.txt"

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  "${PRESPEC}" \
  "${ANALYZER}" \
  "${ROOT}/scripts/score_gse267145_histology_production_v4.py" \
  "${ROOT}/scripts/evaluate_gse267145_histology_predictions.py" \
  "${SCORES}/ARTIFACTS.json" \
  "${OUTCOMES}/ARTIFACTS.json" \
  "${FOLDS}/ARTIFACTS.json" \
  "${SEALED_PRESPEC}/ARTIFACTS.json" \
  slurm/run_gse267145_histology_paired_contrasts_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"

"${PYTHON}" - "${STAGE}" "${PRESPEC_SHA256}" "${SEALED_PRESPEC_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
prespec_sha256, sealed_sha256 = sys.argv[2:]
verify_frozen_tree(root / "contrasts")
receipt = json.loads(
    (root / "contrasts/contrast_receipt.json").read_text(encoding="utf-8")
)
if (
    receipt.get("prespec_contract_sha256") != prespec_sha256
    or receipt.get("sealed_prespec_artifacts_sha256") != sealed_sha256
    or receipt.get("bootstrap_indices_sha256")
    != "fb7ea4b28ef37ee759a70bb3664868f6affa5accc8dadd96edf202238536b2c5"
    or receipt.get("bootstrap_indices_regenerated_and_verified") is not True
    or receipt.get("original_scoring_resamples_reused") is not True
    or receipt.get("frozen_point_estimates_reproduced") != 99
    or receipt.get("declared_contrasts") != 189
    or receipt.get("observed_contrasts") != 169
    or receipt.get("not_applicable_contrasts") != 20
    or receipt.get("per_replicate_difference_vectors_deposited") != 169
    or receipt.get("models_refitted") != 0
    or receipt.get("prediction_bundles_mutated") is not False
    or receipt.get("paired_cluster_bootstrap_used") is not False
    or receipt.get("p_values_calculated") is not False
    or receipt.get("bh_adjustment_calculated") is not False
    or receipt.get("model_ranking_performed") is not False
    or receipt.get("champion_claim_allowed") is not False
    or receipt.get("confirmatory_inference_allowed") is not False
    or receipt.get("between_sex_claim_allowed") is not False
):
    raise SystemExit("paired contrast receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_paired_contrast_campaign",
        "prespec_contract_sha256": prespec_sha256,
        "sealed_prespec_artifacts_sha256": sealed_sha256,
        "bootstrap_indices_sha256": (
            "fb7ea4b28ef37ee759a70bb3664868f6affa5accc8dadd96edf202238536b2c5"
        ),
        "declared_contrasts": 189,
        "observed_contrasts": 169,
        "not_applicable_contrasts": 20,
        "frozen_point_estimates_reproduced": 99,
        "models_refitted": 0,
        "p_values_calculated": False,
        "bh_adjustment_calculated": False,
        "model_ranking_performed": False,
        "status": "passed_paired_development_contrasts",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/contrasts/ARTIFACTS.json"
