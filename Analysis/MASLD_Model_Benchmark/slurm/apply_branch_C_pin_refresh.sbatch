#!/usr/bin/env bash
#SBATCH --job-name=model-check-916
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
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
AUTH="${ROOT}/config/campaigns/authorized_pin_refresh_branch_C_20260825.json"
TARGET="${ROOT}/executions/branch-C-pin-refresh-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# Protocol: fresh bound set immediately before the branch, with its control.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/bound_set.txt"
import json, os, re, subprocess, collections
from pathlib import Path
auth=json.loads(Path("config/campaigns/authorized_pin_refresh_branch_C_20260825.json").read_text())
files=set(auth["branch_closure"])|{p["target_path"] for p in auth["authorized_pins"]}
self_id=os.environ.get("SLURM_JOB_ID","")
ids=[j for j in subprocess.run(["squeue","-u","jameslee","-h","-o","%i"],capture_output=True,text=True).stdout.split()
     if j.split("_")[0]!=self_id]
cfg=collections.Counter(); hits=[]
for jid in ids:
    info=subprocess.run(["scontrol","show","job",jid],capture_output=True,text=True).stdout
    m=re.search(r"Command=(\S+)",info)
    if not m: continue
    c=Path(m.group(1))
    if not c.is_file(): continue
    b=c.read_text(errors="replace")
    for s in re.finditer(r"(scripts/[\w/\-]+\.py)",b):
        q=Path(s.group(1))
        if q.is_file(): b+="\n"+q.read_text(errors="replace")
    for m2 in re.finditer(r"config/[\w/\-]+\.(?:json|toml)",b): cfg[m2.group(0)]+=1
    for f in files:
        if f in b or Path(f).name in b: hits.append((jid,f))
print(f"live jobs {len(ids)}; branch files referenced {len(hits)}")
for h in hits: print("   ",h)
print(f"POSITIVE CONTROL: config files found in live wrappers = {len(cfg)}")
assert len(cfg)>0, "detector did not fire; an empty bound set would mean nothing"
if hits: raise SystemExit("a live job references a branch file; refusing")
PY

PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_refresh_authorized_pins > "${STAGING}/validation/unittest_refresher.txt" 2>&1

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/refresh_authorized_pins.py \
  --root "${ROOT}" --authorization "${AUTH}" \
  --output "${STAGING}/audit/refresh" --apply \
  > "${STAGING}/validation/refresh.stdout.txt" 2>&1

# Branch-complete means: nothing in this branch is left stale, and no pin anywhere
# was stranded by the transaction.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/post_branch_state.txt"
import hashlib, json
from pathlib import Path
prefix="Analysis/MASLD_Model_Benchmark/"
def resolve(i):
    r=i[len(prefix):] if i.startswith(prefix) else i
    c=Path(r); return c if c.is_file() else None
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
                    mism.append(f"{b}: {v}")
        for v in n.values(): walk(v,b)
    elif isinstance(n,list):
        for v in n: walk(v,b)
for b in sorted(Path("config/artifacts/models").glob("*/*.json")):
    walk(json.loads(b.read_text()), f"{b.parent.name}/{b.name}")
for extra in ("config/artifacts/models/cell_bulk_classical_baselines/exposure_audit.json",):
    walk(json.loads(Path(extra).read_text()), extra)
print(f"dead pins across the tree now: {len(mism)} of {checked} checked")
for m in sorted(mism): print("   ",m)
auth=json.loads(Path("config/campaigns/authorized_pin_refresh_branch_C_20260825.json").read_text())
branch={p["bundle_path"] for p in auth["authorized_pins"]}
left=[m for m in mism if any(b.split("/")[-2]+"/"+b.split("/")[-1] in m or b in m for b in branch)]
assert not left, f"branch C left stale pins: {left}"
print("branch C carries no stale pin; nothing in this branch was stranded")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_refresh_authorized_pins tests.unit.test_check_contract_suite_baseline \
  tests.contract.test_registry_files > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_refresh_authorized_pins --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp "${AUTH}" "${STAGING}/audit/"
sha256sum "${AUTH}" "${ROOT}/scripts/refresh_authorized_pins.py" \
  "${ROOT}/config/artifacts/models/cell_bulk_classical_baselines/exposure_audit.json" > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"branch_C_cell_bulk_pin_refresh","branch":"C",
    "pins_refreshed":3,"branch_complete":True,"pins_stranded":0,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
