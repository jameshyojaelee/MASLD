#!/usr/bin/env bash
#SBATCH --job-name=model-check-927
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/register-r10-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# A register about verification must not ship unverified citations.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/register_verification.txt"
import json
from pathlib import Path
UNLOC = "not located during this consolidation"
d = json.loads(Path("config/artifacts/incidents/verification_method_register_20260825.json").read_text())
rules = d["rules"]
ex = [e for r in rules for e in r["worked_examples"]]
loc = [e for e in ex if e["artifact"] != UNLOC]

# 1) every cited output file must resolve on disk
bad = []
for e in loc:
    for path in [p.strip() for p in e["artifact"].split(",")]:
        if not Path(path).exists():
            bad.append(path)
extra = [r["reference_instance_of_the_control"] for r in rules if "reference_instance_of_the_control" in r]
for path in extra:
    if not Path(path).exists(): bad.append(path)
print(f"cited artifacts checked: {len(loc)+len(extra)}, unresolved: {len(bad)}")
print(f"instances still marked as gaps: {len(ex)-len(loc)}")
for b in bad: print("   MISSING", b)
if bad: raise SystemExit("register cites artifacts that do not exist")

# 2) self-counts must be derived, not stale
rc = d["reference_completeness"]
assert rc["rules"] == len(rules), "rule count stale"
assert rc["worked_examples"] == len(ex), "example count stale"
assert rc["instances_with_a_located_artifact"] == len(loc), "located count stale"
assert rc["instances_recorded_without_one"] == len(ex) - len(loc), "unlocated count stale"
print(f"self-counts agree: {len(rules)} rules, {len(ex)} examples, {len(loc)} located, {len(ex)-len(loc)} unlocated")

# 3) every rule must carry the three required parts
missing = [r["id"] for r in rules if not (r.get("rule") and r.get("worked_examples") and r.get("cheap_check"))]
if missing: raise SystemExit("rules missing rule/example/cheap_check: " + ", ".join(missing))
print(f"all {len(rules)} rules carry rule + worked example + cheap check")

# 4) no attribution to individuals. Second person is NOT attribution: an
# instructional register is meant to address its reader ("the field name you
# expect"). An earlier version of this check matched a bare " you " and failed
# on seven legitimate instructional uses -- an over-broad guard, the same class
# of error as a vacuous one, and it fails in the noisier direction.
low = json.dumps(d).lower()
for term in ("team lead", "the user", "i found", "i checked", "agent said",
             "jameslee", "cohortactivation", "as i ", "we decided"):
    assert term not in low, f"register attributes to an individual or session: {term!r}"
print("no attribution to individuals detected (instructional second person is permitted)")
from scripts.check_contract_suite_baseline import PREEXISTING_FAILURES
debt=json.loads(Path("config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json").read_text())
# The debt items array now also holds non-test findings contributed by other
# lanes (a resolved record staleness, a standing recommendation). failing_tests
# counts TESTS, so the cross-check must exclude entries flagged "(not a test)".
openi=[i for i in debt["items"]
       if not i["status"].startswith("RESOLVED")
       and not i["test"].startswith("(not a test)")]
nontest=[i for i in debt["items"] if i["test"].startswith("(not a test)")]
assert len(PREEXISTING_FAILURES)==debt["failing_tests"]==len(openi)==0, "debt and baseline disagree"
assert any("corrections_to_the_brief" == k for k in d), "brief corrections missing from register"
print(f"debt record and baseline agree at {len(PREEXISTING_FAILURES)} tests; "
      f"{len(nontest)} non-test findings carried alongside; brief corrections present")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_check_contract_suite_baseline \
  tests.unit.test_config_write_transaction \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_registry_files | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp config/artifacts/incidents/verification_method_register_20260825.json "${STAGING}/audit/"
cp config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json "${STAGING}/audit/"
sha256sum "${ROOT}/scripts/config_write_transaction.py" \
  "${ROOT}/config/artifacts/incidents/verification_method_register_20260825.json" \
          "${ROOT}/slurm/freeze_verification_method_register.sbatch" > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"verification_method_register","files_modified":0,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
