#!/usr/bin/env bash
#SBATCH --job-name=model-lock-338
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
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
FEATURES=${ROOT}/executions/model-data-332-21100842/view/features
FEATURES_SHA256=1cacb52032beb27969c69ae2d00cfdf3d18b333f1458bb4fb830f5225d79dc6c
BASELINE=${ROOT}/executions/model-cpu-transfer-336-21101362/branch
TRANSCRIPT=${ROOT}/executions/model-cpu-transfer-337-21101402/branch
TARGET=${ROOT}/executions/model-lock-338-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
trap 'status=$?; if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi; exit ${status}' EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${FEATURES_SHA256}" "${FEATURES}/ARTIFACTS.json" | sha256sum --check --strict
BASELINE_SHA256=$(sha256sum "${BASELINE}/ARTIFACTS.json" | cut -d ' ' -f 1)
TRANSCRIPT_SHA256=$(sha256sum "${TRANSCRIPT}/ARTIFACTS.json" | cut -d ' ' -f 1)
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" scripts/lock_gse296875_cell_state_transfer_predictions.py \
  --features "${FEATURES}" --expected-features-sha256 "${FEATURES_SHA256}" \
  --baseline-branch "${BASELINE}" --expected-baseline-sha256 "${BASELINE_SHA256}" \
  --transcriptformer-branch "${TRANSCRIPT}" --expected-transcriptformer-sha256 "${TRANSCRIPT_SHA256}" \
  --output "${STAGE}/lock" \
  > "${STAGE}/validation/lock.stdout.json" 2> "${STAGE}/validation/lock.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum scripts/lock_gse296875_cell_state_transfer_predictions.py \
  slurm/lock_gse296875_cell_state_transfer_predictions_cpu.sbatch \
  "${FEATURES}/ARTIFACTS.json" "${BASELINE}/ARTIFACTS.json" "${TRANSCRIPT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${RUNTIME}/env/bin/python" - <<'PY' "${STAGE}"
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]); meta = verify_frozen_tree(root / "lock")["metadata"]
if meta.get("status") != "locked" or len(meta.get("candidate_roster", [])) != 7:
    raise SystemExit("prediction lock campaign differs")
freeze_tree(root, {
    "artifact_class": "gse296875_cell_state_prediction_lock_campaign",
    "view_id": "gse296875_rna_cell_state_external_development_7500_v1",
    "candidate_roster": meta["candidate_roster"],
    "rows": 7500,
    "evaluator_labels_read": False,
    "evaluator_donor_ids_read": False,
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/lock/ARTIFACTS.json"
