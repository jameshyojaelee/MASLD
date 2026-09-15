#!/usr/bin/env bash
#SBATCH --job-name=model-check-909
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
TARGET="${ROOT}/executions/scgpt-weight-state-correction-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# No live job may be reading the scgpt bundle while it is corrected.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/live_job_recheck.txt"
import os, re, subprocess
from pathlib import Path
# Exclude this job: its own wrapper names scgpt, so it would flag itself.
self_id = os.environ.get("SLURM_JOB_ID", "")
ids = [j for j in subprocess.run(["squeue","-u","jameslee","-h","-o","%i"],
       capture_output=True,text=True).stdout.split() if j.split("_")[0] != self_id]
hits=[]
for jid in ids:
    info = subprocess.run(["scontrol","show","job",jid],capture_output=True,text=True).stdout
    m = re.search(r"Command=(\S+)", info)
    if not m: continue
    cmd = Path(m.group(1))
    if not cmd.is_file(): continue
    blob = cmd.read_text(errors="replace")
    for s in re.finditer(r"(scripts/[\w/\-]+\.py)", blob):
        q=Path(s.group(1))
        if q.is_file(): blob += "\n"+q.read_text(errors="replace")
    if "scgpt" in blob: hits.append(jid)
print(f"live jobs checked: {len(ids)}; referencing scgpt: {len(hits)} {hits}")
if hits: raise SystemExit("a live job references scgpt; refusing to correct")
PY

PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_correct_scgpt_weight_state \
  > "${STAGING}/validation/unittest_corrector.txt" 2>&1

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/correct_scgpt_weight_state.py \
  --root "${ROOT}" --output "${STAGING}/audit/correction" --apply \
  > "${STAGING}/validation/correction.stdout.txt" 2>&1

# The corrected bundle must be true, digests resolved, licence still UNDECLARED,
# and the pinner must rederive against the new bundle.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/post_correction_state.txt"
import hashlib, json
from pathlib import Path
b=json.loads(Path("config/artifacts/models/scgpt/checkpoints.json").read_text())
assert b["weight_content_downloaded"] is True
assert b["weight_license"] == "UNDECLARED", "licence must not be upgraded"
for name,exp in (("scgpt_continual","ad0252a1971e0cd619b7116dbab3177432236c4537225d54280a2aa7e5fe402a"),
                 ("scgpt_whole_human","6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922")):
    a={x["path"]:x for x in b["bundles"][name]["artifacts"]}
    assert a["best_model.pt"]["sha256"]==exp, name
assert "not downloaded" not in b["source_evidence"]
print("bundle corrected; weight_license still UNDECLARED")
live=hashlib.sha256(Path("config/artifacts/models/scgpt/checkpoints.json").read_bytes()).hexdigest()
p=json.loads(Path("config/artifacts/models/scgpt/native_task_disposition_20260825.json").read_text())
pin=[e for e in p["bound_evidence"] if e["path"].endswith("scgpt/checkpoints.json")][0]
assert pin["sha256"]==live, "pinner did not follow the bundle"
assert "weight_content_downloaded is false" not in p["family_specific_blockers"]["detail"]
print("pinner rederives against the corrected bundle; no binding stranded")
for f in ("executions/scgpt-acquisition-21015743/scgpt_continual/best_model.pt",
          "executions/scgpt-acquisition-21015743/scgpt_whole_human/best_model.pt",
          "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/perturbation/checkpoints/scgpt/whole-human/best_model.pt"):
    assert Path(f).is_file(), f"weight body missing - nothing may be deleted: {f}"
print("all 3 weight copies still present; nothing deleted")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_correct_scgpt_weight_state \
  tests.unit.test_check_contract_suite_baseline \
  tests.unit.test_campaign_census_scoping \
  tests.unit.test_refresh_authorized_pins \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_correct_scgpt_weight_state \
  --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

grep -q "test_scgpt_upstream_preflight_and_exposure_are_frozen .* \.\.\. ok" \
  "${STAGING}/validation/unittest.txt" || { echo "SCGPT TEST DID NOT PASS"; exit 1; }
echo "scgpt preflight test now passes"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp config/artifacts/incidents/scgpt_weight_rights_disposition_20260825.json "${STAGING}/audit/"
sha256sum \
  "${ROOT}/config/artifacts/incidents/scgpt_weight_rights_disposition_20260825.json" \
  "${ROOT}/config/artifacts/models/scgpt/checkpoints.json" \
  "${ROOT}/config/artifacts/models/scgpt/native_task_disposition_20260825.json" \
  "${ROOT}/scripts/correct_scgpt_weight_state.py" \
  "${ROOT}/tests/contract/test_registry_files.py" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class": "scgpt_weight_state_correction_and_rights_disposition",
    "weights_deleted": False, "weight_license_changed": False,
    "bundle_corrected": True, "pinner_updated": True,
    "sealed_outcomes_read": False, "models_fit": False, "models_scored": False,
}))
PY
mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
