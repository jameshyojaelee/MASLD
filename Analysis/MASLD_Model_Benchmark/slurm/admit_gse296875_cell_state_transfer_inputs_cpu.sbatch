#!/usr/bin/env bash
#SBATCH --job-name=model-check-336
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
FEATURES=${ROOT}/executions/model-data-332-21100842/view/features
FEATURES_SHA256=1cacb52032beb27969c69ae2d00cfdf3d18b333f1458bb4fb830f5225d79dc6c
VERIFIER=${ROOT}/executions/model-check-331-21100868
VERIFIER_SHA256=9b77301469ae22519b892a1b97a6a392cc9e90e7b857e1c74e2b3541ada45af1
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
BASELINE=${ROOT}/executions/model-cpu-train-503-21097161/predictions
TF_EMBEDDINGS=${ROOT}/executions/model-training-500-21080428
TF_EMBEDDINGS_SHA256=fc4601988491b086d5e89f0ff35f2f73d9c27cf74846a6ff7895efc130c825a5
TF_HEADS=${ROOT}/executions/model-training-21082720
TF_HEADS_SHA256=f0cefd02290c2b2135c3578784076257e2a63213a99e670bf848a47f18a1e27e
TARGET=${ROOT}/executions/model-check-336-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
trap 'status=$?; if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi; exit ${status}' EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
BASELINE_SHA256=$(sha256sum "${BASELINE}/ARTIFACTS.json" | cut -d ' ' -f 1)
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" scripts/admit_gse296875_cell_state_transfer_inputs.py \
  --features "${FEATURES}" --expected-features-sha256 "${FEATURES_SHA256}" \
  --verifier "${VERIFIER}" --expected-verifier-sha256 "${VERIFIER_SHA256}" \
  --source "${SOURCE}" --expected-source-sha256 "${SOURCE_SHA256}" \
  --baseline "${BASELINE}" \
  --tf-embeddings "${TF_EMBEDDINGS}" --expected-tf-embeddings-sha256 "${TF_EMBEDDINGS_SHA256}" \
  --tf-heads "${TF_HEADS}" --expected-tf-heads-sha256 "${TF_HEADS_SHA256}" \
  --output "${STAGE}/authority" \
  > "${STAGE}/validation/admit.stdout.json" 2> "${STAGE}/validation/admit.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum scripts/admit_gse296875_cell_state_transfer_inputs.py \
  slurm/admit_gse296875_cell_state_transfer_inputs_cpu.sbatch \
  "${FEATURES}/ARTIFACTS.json" "${VERIFIER}/ARTIFACTS.json" "${SOURCE}/ARTIFACTS.json" \
  "${BASELINE}/ARTIFACTS.json" "${TF_EMBEDDINGS}/ARTIFACTS.json" "${TF_HEADS}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${RUNTIME}/env/bin/python" - <<'PY' "${STAGE}"
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]); meta = verify_frozen_tree(root / "authority")["metadata"]
if (
    meta.get("baseline_source_fit_census") != 15
    or meta.get("tf_sapiens_source_fit_census_per_head") != 15
    or meta.get("tf_sapiens_source_fit_census_total") != 30
):
    raise SystemExit("transfer input authority differs")
freeze_tree(root, {
    "artifact_class": "gse296875_cell_state_transfer_input_admission_campaign",
    "view_id": "gse296875_rna_cell_state_external_development_7500_v1",
    "baseline_source_fit_census": 15,
    "tf_sapiens_source_fit_census_total": 30,
    "evaluator_labels_read": False,
    "evaluator_donor_ids_read": False,
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/authority/ARTIFACTS.json"
