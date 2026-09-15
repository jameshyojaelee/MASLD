#!/usr/bin/env bash
#SBATCH --job-name=model-check-908
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
TARGET="${ROOT}/executions/job-map-and-debt-revision-${SLURM_JOB_ID}"
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
assert len(PREEXISTING_FAILURES)==debt["failing_tests"]==len(open_items)==6, "debt and baseline disagree"
assert len(resolved)==2
for i in resolved:
    assert i["test"] not in PREEXISTING_FAILURES, f"{i['test']} resolved but still in baseline"
for i in open_items:
    assert i["test"] in PREEXISTING_FAILURES, f"{i['test']} open but absent from baseline"
scg=[i for i in open_items if "scgpt" in i["test"]][0]
assert scg["status"]=="OPEN_HELD_FOR_RIGHTS_DECISION" and "do_not" in scg
print("debt record and baseline constant agree; scgpt held for rights decision")
m=json.loads(Path("config/campaigns/live_job_bundle_bound_set_20260825.json").read_text())
assert m["wrappers_binding_a_closure_file_or_moved_source"]==0
assert m["positive_control"]["config_files_referenced_by_live_wrappers"]>0, "detector never fired"
print(f"job map: {m['live_jobs']} jobs, {m['distinct_wrappers']} wrappers, bound set empty, positive control fired")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_check_contract_suite_baseline \
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
cp config/campaigns/live_job_bundle_bound_set_20260825.json "${STAGING}/audit/"
cp config/campaigns/stale_pin_refresh_plan_20260825.json "${STAGING}/audit/"
cp config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json "${STAGING}/audit/"
sha256sum \
  "${ROOT}/config/campaigns/live_job_bundle_bound_set_20260825.json" \
  "${ROOT}/config/campaigns/stale_pin_refresh_plan_20260825.json" \
  "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class": "job_bundle_bound_set_map_and_debt_revision_3",
    "baseline_failures": 6, "started_at": 8, "cascade_authorized": False,
    "sealed_outcomes_read": False, "models_fit": False, "models_scored": False,
}))
PY
mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
