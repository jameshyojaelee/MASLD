#!/usr/bin/env bash
#SBATCH --job-name=model-check-910
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/debt-revision-4-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/consistency.txt"
import json
from pathlib import Path
from scripts.check_contract_suite_baseline import PREEXISTING_FAILURES
debt=json.loads(Path("config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json").read_text())
open_items=[i for i in debt["items"] if not i["status"].startswith("RESOLVED")]
resolved=[i for i in debt["items"] if i["status"].startswith("RESOLVED")]
print(f"baseline constant: {len(PREEXISTING_FAILURES)}")
print(f"debt record: {debt['failing_tests']} open, {debt['resolved_since_first_record']} resolved")
assert len(PREEXISTING_FAILURES)==debt["failing_tests"]==len(open_items)==5, "debt and baseline disagree"
assert len(resolved)==3
for i in resolved:
    assert i["test"] not in PREEXISTING_FAILURES, f"{i['test']} resolved but still in baseline"
for i in open_items:
    assert i["test"] in PREEXISTING_FAILURES, f"{i['test']} open but absent from baseline"
scg=[i for i in debt["items"] if "scgpt" in i["test"]][0]
assert scg["status"].startswith("RESOLVED"), "scgpt should now be resolved"
assert "geneformer" in json.dumps(debt["open_rights_items"]), "geneformer finding must stay recorded"
b=json.loads(Path("config/artifacts/models/scgpt/checkpoints.json").read_text())
assert b["weight_content_downloaded"] is True and b["weight_license"]=="UNDECLARED"
rd=json.loads(Path("config/artifacts/incidents/scgpt_weight_rights_disposition_20260825.json").read_text())
assert rd["weights_deleted"] is False and rd["standing_constraints"]["champion_ineligible"] is True
assert rd["standing_constraints"]["derivative_weight_redistribution"].startswith("PROHIBITED")
assert rd["exposure_surface"]["total_copies_known"]==3
print("debt record and baseline agree; scgpt corrected; rights disposition intact")

PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_check_contract_suite_baseline \
  tests.unit.test_correct_scgpt_weight_state \
  tests.unit.test_refresh_authorized_pins \
  tests.unit.test_campaign_census_scoping \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp config/artifacts/incidents/scgpt_weight_rights_disposition_20260825.json "${STAGING}/audit/"
cp config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json "${STAGING}/audit/"
sha256sum \
  "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" \
  "${ROOT}/config/artifacts/incidents/scgpt_weight_rights_disposition_20260825.json" \
  "${ROOT}/config/artifacts/models/scgpt/checkpoints.json" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class": "debt_revision_4_after_scgpt_correction",
    "baseline_failures": 5, "started_at": 8, "cascade_authorized": False,
    "geneformer_flagged_unaudited": True,
    "sealed_outcomes_read": False, "models_fit": False, "models_scored": False,
}))
PY
mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
