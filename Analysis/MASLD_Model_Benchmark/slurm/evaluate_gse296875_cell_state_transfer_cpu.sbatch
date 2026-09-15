#!/usr/bin/env bash
#SBATCH --job-name=model-eval-339
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
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
LOCK=${ROOT}/executions/model-lock-338-21101405/lock
EVALUATOR=${ROOT}/executions/model-data-332-21100842/view/evaluator
EVALUATOR_SHA256=d76d571e824e8a76473ce8ecb0ee3279e0b27fff9ca27acb22a4dd34ef2875cc
TARGET=${ROOT}/executions/model-eval-339-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
trap 'status=$?; if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi; exit ${status}' EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EVALUATOR_SHA256}" "${EVALUATOR}/ARTIFACTS.json" | sha256sum --check --strict
LOCK_SHA256=$(sha256sum "${LOCK}/ARTIFACTS.json" | cut -d ' ' -f 1)
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" scripts/evaluate_gse296875_cell_state_transfer.py \
  --prediction-lock "${LOCK}" --expected-prediction-lock-sha256 "${LOCK_SHA256}" \
  --evaluator "${EVALUATOR}" --expected-evaluator-sha256 "${EVALUATOR_SHA256}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluate.stdout.json" 2> "${STAGE}/validation/evaluate.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum scripts/evaluate_gse296875_cell_state_transfer.py \
  slurm/evaluate_gse296875_cell_state_transfer_cpu.sbatch \
  "${LOCK}/ARTIFACTS.json" "${EVALUATOR}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${RUNTIME}/env/bin/python" - <<'PY' "${STAGE}"
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]); meta = verify_frozen_tree(root / "evaluation")["metadata"]
if (
    meta.get("biological_replication_unit") != "donor"
    or meta.get("donors") != 39
    or len(meta.get("candidate_roster", [])) != 7
    or meta.get("external_or_champion_claim_allowed") is not False
):
    raise SystemExit("transfer evaluator campaign differs")
freeze_tree(root, {
    "artifact_class": "gse296875_cell_state_transfer_evaluation_campaign",
    "view_id": "gse296875_rna_cell_state_external_development_7500_v1",
    "candidate_roster": meta["candidate_roster"],
    "rows": 7500,
    "donors": 39,
    "biological_replication_unit": "donor",
    "project_exposed_development_only": True,
    "external_or_champion_claim_allowed": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/evaluation/ARTIFACTS.json"
