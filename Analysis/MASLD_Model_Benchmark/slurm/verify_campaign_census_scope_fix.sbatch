#!/usr/bin/env bash
#SBATCH --job-name=model-check-905
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
TARGET="${ROOT}/executions/campaign-census-scope-fix-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# Re-derive the campaign partition the fix relies on.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/campaign_partition.txt"
import tomllib
from pathlib import Path
from tests.contract.test_registry_files import (
    CAMPAIGN_MODEL_DISPATCH_KEYS, campaign_dispatch_keys, references_registry_model_id)
bound=inert=partial=0
for p in sorted(Path("config/campaigns").glob("*.toml")):
    d=tomllib.load(open(p,'rb'))
    k=campaign_dispatch_keys(d)
    if not k:
        assert not references_registry_model_id(d), f"{p.name} names a registry model"
        inert+=1
    elif k==set(CAMPAIGN_MODEL_DISPATCH_KEYS): bound+=1
    else: partial+=1; print("PARTIAL:",p.name,sorted(set(CAMPAIGN_MODEL_DISPATCH_KEYS)-k))
print(f"bound={bound} inert={inert} partial={partial}")
assert (bound,inert,partial)==(18,5,0), "campaign partition changed"
print("campaign partition verified: 18 bind, 5 inert, 0 partial")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_campaign_census_scoping \
  tests.unit.test_check_contract_suite_baseline \
  tests.unit.test_build_pxd051911_activation_readiness \
  tests.unit.test_build_gse256398_terminal_disposition \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_campaign_census_scoping \
  --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

# The fixed test must now be reported as passing, not merely absent.
grep -q "test_every_campaign_binds_the_frozen_model_census .* \.\.\. ok" \
  "${STAGING}/validation/unittest.txt" \
  || { echo "FIXED TEST DID NOT REPORT ok"; exit 1; }
echo "census-binding test now passes"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
sha256sum \
  "${ROOT}/tests/contract/test_registry_files.py" \
  "${ROOT}/tests/unit/test_campaign_census_scoping.py" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" \
  "${ROOT}/slurm/verify_campaign_census_scope_fix.sbatch" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class": "campaign_census_scope_fix",
    "resolved_contract_failure": "test_every_campaign_binds_the_frozen_model_census",
    "baseline_failures_before": 8, "baseline_failures_after": 7,
    "campaign_data_edited": False,
    "sealed_outcomes_read": False, "models_fit": False, "models_scored": False,
}))
PY
mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
