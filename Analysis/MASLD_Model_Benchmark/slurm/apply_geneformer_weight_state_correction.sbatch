#!/usr/bin/env bash
#SBATCH --job-name=model-check-912
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
TARGET="${ROOT}/executions/geneformer-weight-state-correction-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/live_job_recheck.txt"
import os, re, subprocess
from pathlib import Path
self_id = os.environ.get("SLURM_JOB_ID","")
ids=[j for j in subprocess.run(["squeue","-u","jameslee","-h","-o","%i"],capture_output=True,text=True).stdout.split()
     if j.split("_")[0]!=self_id]
hits=[]
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
    if "geneformer" in b: hits.append(jid)
print(f"live jobs checked: {len(ids)}; referencing geneformer: {len(hits)} {hits}")
if hits: raise SystemExit("a live job references geneformer; refusing to correct")
PY

PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_correct_geneformer_weight_state \
  > "${STAGING}/validation/unittest_corrector.txt" 2>&1

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/correct_geneformer_weight_state.py \
  --root "${ROOT}" --output "${STAGING}/audit/correction" --apply \
  > "${STAGING}/validation/correction.stdout.txt" 2>&1

PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/post_correction_state.txt"
import hashlib, json
from pathlib import Path
b=json.loads(Path("config/artifacts/models/geneformer/checkpoints.json").read_text())
assert b["weight_content_downloaded"] is True
assert b["weight_license"]=="Apache-2.0", "licence must not be downgraded"
assert b["code_license"]=="Apache-2.0"
assert "license" not in b, "singular key should be normalised away"
print("bundle corrected; weight_license = Apache-2.0, NOT downgraded")
live=hashlib.sha256(Path("config/artifacts/models/geneformer/checkpoints.json").read_bytes()).hexdigest()
p=json.loads(Path("config/artifacts/models/geneformer/native_task_disposition_20260825.json").read_text())
pin=[e for e in p["bound_evidence"] if e["path"].endswith("geneformer/checkpoints.json")][0]
assert pin["sha256"]==live, "pinner did not follow the bundle"
ri=p["recorded_inconsistencies"][0]
assert ri["action"]=="corrected_2026_08_25" and ri["recorded_value"] is False
print("pinner rederives; deferred item marked closed with its original observation intact")
# every weight body must still exist, including the abandoned staging copies
for d in ("executions/geneformer-acquisition-21014827/source",
          "executions/geneformer-acquisition-21014704.staging/source",
          "executions/geneformer-acquisition-21014802.staging/source"):
    for n in ("Geneformer-V1-10M/model.safetensors","Geneformer-V2-104M/model.safetensors",
              "Geneformer-V2-316M/model.safetensors"):
        assert (Path(d)/n).is_file(), f"weight body missing - nothing may be deleted: {d}/{n}"
print("all 9 weight body copies still present; nothing deleted")
# the withdrawn finding must remain visible, not silently removed
for f in ("config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json",
          "config/artifacts/incidents/scgpt_weight_rights_disposition_20260825.json"):
    t=Path(f).read_text()
    assert "WITHDRAWN" in t and "what_was_claimed" in t, f"withdrawn finding not visible in {f}"
print("withdrawn finding remains visible in both records")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_correct_geneformer_weight_state \
  tests.unit.test_correct_scgpt_weight_state \
  tests.unit.test_check_contract_suite_baseline \
  tests.unit.test_campaign_census_scoping \
  tests.unit.test_refresh_authorized_pins \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_correct_geneformer_weight_state \
  --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp config/artifacts/incidents/geneformer_weight_state_correction_20260825.json "${STAGING}/audit/"
sha256sum \
  "${ROOT}/config/artifacts/incidents/geneformer_weight_state_correction_20260825.json" \
  "${ROOT}/config/artifacts/models/geneformer/checkpoints.json" \
  "${ROOT}/config/artifacts/models/geneformer/native_task_disposition_20260825.json" \
  "${ROOT}/config/artifacts/incidents/contract_suite_preexisting_debt_20260825.json" \
  "${ROOT}/scripts/correct_geneformer_weight_state.py" \
  > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"geneformer_weight_state_correction",
    "rights_gap_found":False,"weight_license":"Apache-2.0","weight_license_downgraded":False,
    "weights_deleted":False,"digests_corrected":0,"earlier_finding_withdrawn":True,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
