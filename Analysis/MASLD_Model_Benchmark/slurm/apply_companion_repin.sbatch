#!/usr/bin/env bash
#SBATCH --job-name=model-check-921
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
AUTH="${ROOT}/config/campaigns/authorized_pin_refresh_companion_repin_20260825.json"
TARGET="${ROOT}/executions/companion-repin-${SLURM_JOB_ID}"; STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/refresh_authorized_pins.py \
  --root "${ROOT}" --authorization "${AUTH}" --output "${STAGING}/audit/repin" --apply \
  > "${STAGING}/validation/repin.stdout.txt" 2>&1

# Whole-tree pin check by digest value, not by shape.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/whole_tree_pins.txt"
import hashlib, json
from pathlib import Path
prefix="Analysis/MASLD_Model_Benchmark/"
checked=0; mism=[]
def resolve(v):
    r=v[len(prefix):] if v.startswith(prefix) else v
    p=Path(r); return p if p.is_file() else None
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
                if hashlib.sha256(s.read_bytes()).hexdigest()!=n[dk]: mism.append(f"{b}: {v}")
        for v in n.values(): walk(v,b)
    elif isinstance(n,list):
        for v in n: walk(v,b)
for b in sorted(Path("config/artifacts/models").glob("*/*.json")):
    walk(json.loads(b.read_text()), f"{b.parent.name}/{b.name}")
walk(json.loads(Path("config/typed_evidence_graph_tournament_preflight.json").read_text()),
     "typed_evidence_graph_tournament_preflight.json")
print(f"dead pins across the tree: {len(mism)} of {checked} checked")
for m in sorted(mism): print("   ",m)
assert not mism, "pins still dead"
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_refresh_pin_companion_metadata tests.unit.test_refresh_authorized_pins \
  tests.unit.test_check_contract_suite_baseline tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_registry_files | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp "${AUTH}" "${STAGING}/audit/"
sha256sum "${AUTH}" "${ROOT}/scripts/refresh_pin_companion_metadata.py" > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"companion_rewrite_repin","pins_repinned":6,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
