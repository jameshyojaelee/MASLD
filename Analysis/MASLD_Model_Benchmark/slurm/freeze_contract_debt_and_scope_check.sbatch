#!/usr/bin/env bash
#SBATCH --job-name=model-check-904
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
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
TARGET="${ROOT}/executions/contract-debt-and-scope-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# Re-derive every count the debt record asserts, so the record cannot rot quietly.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/debt_rederivation.txt"
import hashlib, json, tomllib, collections
from pathlib import Path
ROOT = Path("."); CONFIG = ROOT / "config"
rec = json.loads((CONFIG/"artifacts/incidents/contract_suite_preexisting_debt_20260825.json").read_text())
problems = []

# 1) stale repo-path pins, using the contract test's exact resolve and key order
prefix = "Analysis/MASLD_Model_Benchmark/"
checked = 0; mism = []
def resolve(i):
    r = i[len(prefix):] if i.startswith(prefix) else i
    c = ROOT / r
    return c if c.is_file() else None
def walk(n, b):
    global checked
    if isinstance(n, dict):
        dk = next((k for k in ("sha256","authority_sha256")
                   if isinstance(n.get(k), str) and len(n[k]) == 64), None)
        if dk:
            for pk in ("identity","authority_path","path"):
                v = n.get(pk)
                if not isinstance(v, str): continue
                s = resolve(v)
                if s is None: continue
                checked += 1
                if hashlib.sha256(s.read_bytes()).hexdigest() != n[dk]:
                    mism.append(v)
        for v in n.values(): walk(v, b)
    elif isinstance(n, list):
        for v in n: walk(v, b)
for b in sorted((CONFIG/"artifacts/models").glob("*/*.json")):
    walk(json.loads(b.read_text(encoding="utf-8")), b)
print(f"stale repo-path pins: {len(mism)} of {checked} checked")
if len(mism) != 16: problems.append(f"stale pin count {len(mism)} != recorded 16")
for t, c in collections.Counter(mism).most_common(): print(f"  {c}x {t}")

# 2) campaigns lacking [selection]
missing = [p.name for p in sorted((CONFIG/"campaigns").glob("*.toml"))
           if "selection" not in tomllib.load(open(p,'rb'))]
print(f"toml campaigns lacking [selection]: {len(missing)}")
for m in missing: print(f"  {m}")
if len(missing) != 5: problems.append(f"campaign count {len(missing)} != recorded 5")

# 3) The census item is quoted from the contract test, not re-derived here. A
# plain set difference does NOT reproduce the test's construction: group_bindings
# holds 17 group ids under the key bundle_id, not per-model bundle ids. Report
# the raw counts only, and assert nothing from them.
cens = tomllib.load(open(CONFIG/"evaluation/model_audit_bindings.toml",'rb'))
on_disk = {p.name for p in (CONFIG/"artifacts/models").iterdir() if p.is_dir()}
print(f"bundle dirs on disk: {len(on_disk)}; "
      f"expected_artifact_bundle_count: {cens.get('expected_artifact_bundle_count')}; "
      f"group_bindings entries: {len(cens.get('group_bindings', []))}")
print("census orphan item: quoted from the contract test, not re-derived here")

# 4) nothing this campaign touched appears in the stale set
for guarded in ("config/datasets/pxd051911.toml", "config/datasets/gse256398.toml",
                "config/evaluation/family_native_tournament.toml"):
    if any(guarded in m for m in mism):
        problems.append(f"campaign-touched file in stale set: {guarded}")
print("campaign-touched files in stale set: none" if not problems else "")

if problems:
    raise SystemExit("DEBT RECORD OUT OF DATE:\n" + "\n".join(problems))
print("debt record rederives")
PY

# Read-only Resource-scope confirmation. Records the result; never acts on it.
set +e
python "${REPO}/scripts/manuscript/validate_resource_scope.py" \
  > "${STAGING}/validation/validate_resource_scope.stdout.txt" \
  2> "${STAGING}/validation/validate_resource_scope.stderr.txt"
echo $? > "${STAGING}/validation/validate_resource_scope.exit"
set -e
echo "validate_resource_scope exit: $(cat "${STAGING}/validation/validate_resource_scope.exit")"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" "${STAGING}/audit/"
sha256sum \
  "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  "${REPO}/scripts/manuscript/validate_resource_scope.py" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class": "contract_suite_debt_and_resource_scope_check",
    "failing_contract_tests": 8,
    "resource_scope_check": "read_only",
    "sealed_outcomes_read": False,
    "models_fit": False,
    "models_scored": False,
}))
PY
mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
