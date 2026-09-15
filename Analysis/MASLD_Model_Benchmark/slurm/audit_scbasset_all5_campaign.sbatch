#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-033
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SCRIPT=${ROOT}/scripts/audit_scbasset_all5_campaign.py
TEST=${ROOT}/tests/unit/test_audit_scbasset_all5_campaign.py
AGGREGATE=${ROOT}/executions/scbasset-all5-valid-meta-aggregate-21070675
AGGREGATE_SHA256=79ecbd2d473c84ebb75e5e5476dcb4059ac27fe1d061263c733e9d35dc1f9a3f
TARGET=${ROOT}/executions/model-cpu-audit-033-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite scBasset campaign audit.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' "${AGGREGATE_SHA256}" "${AGGREGATE}/ARTIFACTS.json" \
  | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - <<'PY'
from pathlib import Path
from masld_bench.artifacts import verify_frozen_tree

root = Path.cwd()
paths = (
    "executions/scbasset-diagonal-matrix-gate-21067170",
    "executions/scbasset-crossfold-evaluator-source-21070148",
    "executions/scbasset-crossfold-evaluator-source-21070472",
    "executions/scbasset-fold1to4-chain-21069969",
    "executions/scbasset-fold1to4-prediction-lock-21070492",
    "executions/scbasset-donor0-genomic0-seed11-valid-evaluation-21065092",
    "executions/scbasset-fold1to4-valid-crossfold-evaluation-21070600",
    "executions/scbasset-all5-valid-meta-aggregate-21070675",
)
for value in paths:
    verify_frozen_tree(root / value)
for fold, suffix in ((0, "21064455"), (1, "21067171"), (2, "21067172"), (3, "21067173"), (4, "21067174")):
    separator = "-" if fold == 0 else "_"
    verify_frozen_tree(root / f"executions/scbasset-donor{fold}{separator}genomic{fold}-inputs-{suffix}")
models = (
    "executions/scbasset-donor0-genomic0-seed11-21064694",
    "executions/scbasset-donor1_genomic1-seed20260824-21069969",
    "executions/scbasset-donor2_genomic2-seed20260824-21069969",
    "executions/scbasset-donor3_genomic3-seed20260824-21069969",
    "executions/scbasset-donor4_genomic4-seed20260824-21069969",
)
predictions = (
    "executions/scbasset-donor0-genomic0-seed11-valid-21064859",
    "executions/scbasset-donor1_genomic1-seed20260824-valid-21069969",
    "executions/scbasset-donor2_genomic2-seed20260824-valid-21069969",
    "executions/scbasset-donor3_genomic3-seed20260824-valid-21069969",
    "executions/scbasset-donor4_genomic4-seed20260824-valid-21069969",
)
for value in models + predictions:
    verify_frozen_tree(root / value)
PY
PYTHONPATH=${ROOT}:${ROOT}/src python -m py_compile "${SCRIPT}" "${TEST}"
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_audit_scbasset_all5_campaign \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH=${ROOT}:${ROOT}/src python "${SCRIPT}" \
  --root "${ROOT}" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json"
sha256sum \
  "${SCRIPT}" "${TEST}" "$0" \
  executions/scbasset-crossfold-evaluator-source-21070148/ARTIFACTS.json \
  executions/scbasset-crossfold-evaluator-source-21070472/ARTIFACTS.json \
  executions/scbasset-fold1to4-prediction-lock-21070492/ARTIFACTS.json \
  executions/scbasset-all5-valid-meta-aggregate-21070675/ARTIFACTS.json \
  > "${STAGE}/source.sha256"
PYTHONPATH=${ROOT}/src python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/campaign_audit.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass"
    or receipt.get("outer_folds_present") != [0, 1, 2, 3, 4]
    or receipt.get("donor_safe_five_fold_artifact_campaign_complete") is not True
    or receipt.get("uniform_prediction_lock_five_fold_campaign_complete") is not False
    or receipt.get("fixed_seed_three_seed_screen_complete") is not False
    or receipt.get("development_gate_passed") is not False
    or receipt.get("gpu_production_authorized") is not False
    or receipt.get("raw_outcome_authority_opened") is not False
    or receipt.get("outcomes_used_for_fitting") is not False
    or receipt.get("sealed_data_read") is not False
):
    raise SystemExit("scBasset campaign audit contract differs")
freeze_tree(root, {
    "artifact_class": "scbasset_all5_campaign_audit",
    "model_id": "scbasset",
    "dataset_id": "gse296875",
    "outer_folds": [0, 1, 2, 3, 4],
    "donor_safe_artifact_campaign_complete": True,
    "fixed_seed_three_seed_screen_complete": False,
    "development_gate_passed": False,
    "gpu_production_authorized": False,
    "sealed_data_read": False,
    "status": "pass",
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
