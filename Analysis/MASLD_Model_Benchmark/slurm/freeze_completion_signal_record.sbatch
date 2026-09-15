#!/usr/bin/env bash
#SBATCH --job-name=model-check-913
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
TARGET="${ROOT}/executions/freeze-completion-signal-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# Re-derive the record's own evidence so it cannot be frozen while stale.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/signal_rederivation.txt"
import json
from pathlib import Path
r = json.loads(Path("config/artifacts/incidents/freeze_completion_signal_definition_20260825.json").read_text())
ev = r["evidence"]["verified_directly"]
problems = []
for e in ev:
    c = Path(e["candidate_present"].split(",")[0].strip())
    print(f"{e['version']}: candidate {'present' if c.is_dir() else 'MISSING'}  {c}")
    if not c.is_dir(): problems.append(f"{e['version']} candidate missing")
    for line in e.get("subsequent_job_outcomes", []):
        p = Path(line.split()[0])
        has = (p / "COMPLETE").is_file()
        print(f"    {p.name}: COMPLETE={'yes' if has else 'no'}")
        if has: problems.append(f"{p} unexpectedly has COMPLETE")
    if "also_passed" in e:
        p = Path(e["also_passed"].split()[0])
        print(f"    {p.name}: COMPLETE={'yes' if (p/'COMPLETE').is_file() else 'no'} (expected yes)")
        if not (p/"COMPLETE").is_file(): problems.append(f"{p} lost its COMPLETE")
# Every entry must declare a verification status; the unverified slot is kept
# deliberately so a future addition cannot inherit the credibility around it.
u = r["evidence"]["asserted_but_not_verified_here"]
assert u["retained_deliberately"] is True, "the unverified slot must be retained"
assert u["entries"] == [], "an unverified entry appeared; it must be re-checked or labelled"
assert r["evidence"]["verified_case_count"] == len(ev), "verified count disagrees with the entries"
print(f"\n{len(ev)} verified cases; unverified slot retained and empty")
if problems: raise SystemExit("record is stale:\n" + "\n".join(problems))
print("freeze-completion signal record re-derives")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_check_contract_suite_baseline \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_registry_files | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp config/artifacts/incidents/freeze_completion_signal_definition_20260825.json "${STAGING}/audit/"
sha256sum "${ROOT}/config/artifacts/incidents/freeze_completion_signal_definition_20260825.json" \
          "${ROOT}/slurm/freeze_completion_signal_record.sbatch" > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"freeze_completion_signal_definition",
    "files_modified":0,"cascade_started":False,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
