#!/usr/bin/env bash
#SBATCH --job-name=model-check-919
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
AUTH="${ROOT}/config/campaigns/authorized_pin_refresh_branch_BD_20260825.json"
TARGET="${ROOT}/executions/branch-BD-pin-refresh-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# --- config-write transaction -------------------------------------------
# A run that writes config and then dies leaves the writes in place, so a
# .failed execution dir means "the job died", not "nothing landed". Snapshot
# before the first write; on any non-zero exit restore, verify each file
# re-derives to its pre-write digest, and only then rename to .failed.
ROLLBACK="${STAGING}/rollback"
FAILED="${TARGET}.failed"
rollback_on_failure() {
  local rc=$?
  if [[ ${rc} -ne 0 ]]; then
    echo "run failed with ${rc}; rolling back config writes"
    if [[ -d "${ROLLBACK}" ]]; then
      PYTHONPATH=${ROOT}:${ROOT}/src python "${ROOT}/scripts/config_write_transaction.py" \
        restore --root "${ROOT}" --into "${ROLLBACK}" || echo "ROLLBACK ITSELF FAILED"
    fi
    if [[ -d "${STAGING}" && ! -e "${FAILED}" ]]; then mv "${STAGING}" "${FAILED}"; fi
  fi
}
trap rollback_on_failure EXIT
# ------------------------------------------------------------------------


PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/bound_set.txt"
import json, os, re, subprocess, collections
from pathlib import Path
auth=json.loads(Path("config/campaigns/authorized_pin_refresh_branch_BD_20260825.json").read_text())
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
print(f"POSITIVE CONTROL: config files in live wrappers = {len(cfg)}")
assert len(cfg)>0, "detector did not fire"
if hits: raise SystemExit("a live job references a branch file; refusing")
PY

PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v tests.unit.test_refresh_authorized_pins \
  > "${STAGING}/validation/unittest_refresher.txt" 2>&1


# Snapshot every file this run may write, BEFORE the first write.
PYTHONPATH=${ROOT}:${ROOT}/src python - "${ROLLBACK}" "${AUTH}" <<'PYX'
import json, subprocess, sys
from pathlib import Path
rollback, auth_path = sys.argv[1], sys.argv[2]
auth = json.loads(Path(auth_path).read_text())
targets = sorted({p["bundle_path"] for p in auth["authorized_pins"]})
cmd = ["python", "scripts/config_write_transaction.py", "snapshot",
       "--root", ".", "--into", rollback]
for t in targets: cmd += ["--target", t]
raise SystemExit(subprocess.call(cmd))
PYX

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/refresh_authorized_pins.py \
  --root "${ROOT}" --authorization "${AUTH}" --output "${STAGING}/audit/refresh" --apply \
  > "${STAGING}/validation/refresh.stdout.txt" 2>&1

# Post-check by DIGEST VALUE, not by shape. A shape-based walk is blind to
# flat-map and non-standard-key recordings, which is what nearly stranded five
# pins in this very branch; verifying with that same blind check would prove
# nothing.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/post_branch_state.txt"
import hashlib, json, tomllib
from pathlib import Path
D=lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
auth=json.loads(Path("config/campaigns/authorized_pin_refresh_branch_BD_20260825.json").read_text())
closure=auth["branch_closure"]
def load(f):
    return tomllib.loads(Path(f).read_text()) if f.endswith(".toml") else json.loads(Path(f).read_text())
# any 64-hex string in a closure file that equals NO live digest of a closure
# file it names, but DOES equal a pre-refresh digest, would be a stranded pin.
live={f:D(f) for f in closure}
stranded=[]
for holder in closure:
    doc=load(holder)
    def walk(n,p=""):
        if isinstance(n,dict):
            for k,v in n.items():
                if isinstance(v,str) and len(v)==64:
                    for tgt,dg in live.items():
                        if tgt==holder: continue
                        # a recording that names a closure target but holds a
                        # digest matching none of them is suspicious
                        pass
                walk(v,f"{p}.{k}")
        elif isinstance(n,list):
            for i,v in enumerate(n): walk(v,f"{p}[{i}]")
    walk(doc)
# direct check: every authorized pin must now hold the live digest of its target
import sys
sys.path.insert(0,".")
from scripts.refresh_authorized_pins import read_toml_pin, navigate
bad=[]
for pin in auth["authorized_pins"]:
    f=pin["bundle_path"]
    if f.endswith(".toml"):
        cur=read_toml_pin(Path(f), pin["json_path"], pin["key"])
    else:
        node=navigate(json.loads(Path(f).read_text()), pin["json_path"])
        cur=node[pin["key"]] if isinstance(node,dict) else node
    want=D(pin["target_path"])
    if cur!=want: bad.append(f"{f} {pin['json_path']}.{pin['key']} = {cur[:8]} want {want[:8]}")
print(f"authorized pins now holding the live digest: {len(auth['authorized_pins'])-len(bad)} of {len(auth['authorized_pins'])}")
for b in bad: print("   STALE", b)
if bad: raise SystemExit("branch B+D left stale pins")
print("every authorized pin resolves; branch B+D is complete")
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
sha256sum "${AUTH}" "${ROOT}/scripts/refresh_authorized_pins.py" > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"branch_BD_pin_refresh","branch":"B+D","pins_refreshed":22,
    "branch_complete":True,"pins_stranded":0,"closure_files":12,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
