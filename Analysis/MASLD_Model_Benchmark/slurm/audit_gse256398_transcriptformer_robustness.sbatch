#!/usr/bin/env bash
#SBATCH --job-name=model-audit-256
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
: "${EMBEDDINGS_ROOT:?EMBEDDINGS_ROOT is required}"
: "${EMBEDDINGS_SHA256:?EMBEDDINGS_SHA256 is required}"
: "${PREDICTION_LOCK_ROOT:?PREDICTION_LOCK_ROOT is required}"
: "${PREDICTION_LOCK_SHA256:?PREDICTION_LOCK_SHA256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
METADATA=${ROOT}/executions/gse256398-donor-metadata-21069678
TARGET=${ROOT}/executions/gse256398-transcriptformer-descriptive-audit-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 TranscriptFormer descriptive audit.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_audit_gse256398_transcriptformer_robustness \
  > "${STAGE}/validation/unittest.txt" 2>&1
"${RUNTIME}/env/bin/python" scripts/audit_gse256398_transcriptformer_robustness.py \
  --embeddings-root "${EMBEDDINGS_ROOT}" \
  --embeddings-sha256 "${EMBEDDINGS_SHA256}" \
  --prediction-lock-root "${PREDICTION_LOCK_ROOT}" \
  --prediction-lock-sha256 "${PREDICTION_LOCK_SHA256}" \
  --metadata-root "${METADATA}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.txt" \
  2> "${STAGE}/validation/audit.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  scripts/audit_gse256398_transcriptformer_robustness.py \
  tests/unit/test_audit_gse256398_transcriptformer_robustness.py \
  slurm/audit_gse256398_transcriptformer_robustness.sbatch \
  "${EMBEDDINGS_ROOT}/ARTIFACTS.json" \
  "${PREDICTION_LOCK_ROOT}/ARTIFACTS.json" \
  "${METADATA}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${RUNTIME}/env/bin/python" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_descriptive_only_no_executable_task_spec"
    or receipt["unit_of_replication"] != "donor"
    or receipt["masld_relevant_donors"] != 17
    or receipt["alcohol_associated_ood_donors"] != 9
    or not receipt["prediction_locked_before_metadata_join"]
    or receipt["sealed_outcomes_read"]
    or receipt["barcode_cell_labels_read"]
    or receipt["supervised_accuracy_claim_allowed"]
    or receipt["model_selection_allowed"]
    or receipt["alcohol_donors_used_as_masld_negative_controls"]
    or receipt["confirmatory_transfer_claim_allowed"]
):
    raise SystemExit("GSE256398 descriptive-only boundary differs")
freeze_tree(root, {
    "artifact_class": "gse256398_transcriptformer_descriptive_robustness_audit",
    "dataset_id": "gse256398",
    "model_id": "transcriptformer_tf_sapiens",
    "unit_of_replication": "donor",
    "supervised_accuracy_claim_allowed": False,
    "model_selection_allowed": False,
    "sealed_outcomes_read": False,
    "status": "passed_descriptive_only",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
