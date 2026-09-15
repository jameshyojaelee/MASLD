#!/usr/bin/env bash
#SBATCH --job-name=model-check-907
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
TARGET="${ROOT}/executions/isolated-pin-refresh-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# Re-confirm, at execution time, that no live job touches the three bundles.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/live_job_recheck.txt"
import json, re, subprocess
from pathlib import Path
auth = json.loads(Path("config/campaigns/authorized_pin_refresh_isolated_20260825.json").read_text())
names = {Path(p["bundle_path"]).parent.name for p in auth["authorized_pins"]}
ids = subprocess.run(["squeue","-u","jameslee","-h","-o","%i"],capture_output=True,text=True).stdout.split()
hits = []
for jid in ids:
    info = subprocess.run(["scontrol","show","job",jid],capture_output=True,text=True).stdout
    m = re.search(r"Command=(\S+)", info)
    if not m: continue
    cmd = Path(m.group(1))
    if not cmd.is_file(): continue
    blob = cmd.read_text(errors="replace")
    for s in re.finditer(r"(scripts/[\w/\-]+\.py)", blob):
        q = Path(s.group(1))
        if q.is_file(): blob += "\n" + q.read_text(errors="replace")
    for n in names:
        if n in blob: hits.append((jid, n))
print(f"live jobs checked: {len(ids)}; references to the authorized bundles: {len(hits)}")
for h in hits: print("  ", h)
if hits: raise SystemExit("a live job references an authorized bundle; refusing to refresh")
print("no live job references the authorized bundles")
PY

PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_refresh_authorized_pins \
  > "${STAGING}/validation/unittest_refresher.txt" 2>&1

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/refresh_authorized_pins.py \
  --root "${ROOT}" \
  --authorization "${ROOT}/config/campaigns/authorized_pin_refresh_isolated_20260825.json" \
  --output "${STAGING}/audit/refresh" \
  --apply \
  > "${STAGING}/validation/refresh.stdout.txt" 2>&1

# The three pins must now rederive, and the cascade must be untouched.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/post_refresh_pin_state.txt"
import hashlib, json
from pathlib import Path
ROOT=Path("."); prefix="Analysis/MASLD_Model_Benchmark/"
def resolve(i):
    r=i[len(prefix):] if i.startswith(prefix) else i
    c=ROOT/r
    return c if c.is_file() else None
checked=0; mism=[]
def walk(n,b):
    global checked
    if isinstance(n,dict):
        dk=next((k for k in ("sha256","authority_sha256") if isinstance(n.get(k),str) and len(n[k])==64),None)
        if dk:
            for pk in ("identity","authority_path","path"):
                v=n.get(pk)
                if not isinstance(v,str): continue
                s=resolve(v)
                if s is None: continue
                checked+=1
                if hashlib.sha256(s.read_bytes()).hexdigest()!=n[dk]:
                    mism.append(f"{b.parent.name}/{b.name}: {v}")
        for v in n.values(): walk(v,b)
    elif isinstance(n,list):
        for v in n: walk(v,b)
for b in sorted((ROOT/"config/artifacts/models").glob("*/*.json")):
    walk(json.loads(b.read_text(encoding="utf-8")),b)
print(f"dead pins now: {len(mism)} of {checked} checked (was 16)")
for m in sorted(mism): print("  ", m)
auth=json.loads(Path("config/campaigns/authorized_pin_refresh_isolated_20260825.json").read_text())
for p in auth["authorized_pins"]:
    tag=f"{Path(p['bundle_path']).parent.name}/{Path(p['bundle_path']).name}"
    assert not any(m.startswith(tag+":") for m in mism), f"{tag} still stale"
assert len(mism)==13, f"expected 13 remaining, got {len(mism)}"
print("the 3 authorized pins now rederive; the other 13 remain, cascade untouched")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_refresh_authorized_pins \
  tests.unit.test_campaign_census_scoping \
  tests.unit.test_check_contract_suite_baseline \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_refresh_authorized_pins \
  --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

grep -q "test_uce_upstream_preflight_and_exposure_are_frozen .* \.\.\. ok" \
  "${STAGING}/validation/unittest.txt" || { echo "UCE TEST DID NOT PASS"; exit 1; }
echo "uce preflight test now passes"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
sha256sum \
  "${ROOT}/config/campaigns/authorized_pin_refresh_isolated_20260825.json" \
  "${ROOT}/scripts/refresh_authorized_pins.py" \
  "${ROOT}/tests/unit/test_refresh_authorized_pins.py" \
  "${ROOT}/tests/contract/test_registry_files.py" \
  "${ROOT}/slurm/apply_isolated_pin_refresh.sbatch" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class": "isolated_pin_refresh_and_uce_test_correction",
    "pins_refreshed": 3, "pins_remaining_dead": 13, "cascade_started": False,
    "scgpt_bundle_touched": False,
    "sealed_outcomes_read": False, "models_fit": False, "models_scored": False,
}))
PY
mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
