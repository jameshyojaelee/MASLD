#!/usr/bin/env bash
#SBATCH --job-name=model-check-345
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
TARGET="${ROOT}/executions/model-check-345-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# The record claims the two pinned sites are UNCHANGED and that no pin moved.
# Re-derive that claim so the record cannot be frozen while stale.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/pin_and_site_rederivation.txt"
import hashlib, json
from pathlib import Path

ck = Path("config/artifacts/models/mpralegnet/checkpoints.json")
live = hashlib.sha256(ck.read_bytes()).hexdigest()
expected = "12f9418d229274b8f14feceed1b1bf4aa918a040326ce85193ad81bf0073880e"
print(f"mpralegnet/checkpoints.json live sha256 = {live}")
problems = []
if live != expected:
    problems.append("checkpoints.json changed; sites A/C may have landed -- re-derive the pins")

# Both pins must still be alive, i.e. still equal to the live digest.
pins = [
    Path("config/gse281364_comparator_reconciliation.toml"),
    Path("config/local_variant_family_reconciliation.json"),
]
for p in pins:
    text = p.read_text(encoding="utf-8")
    alive = live in text
    print(f"  pin in {p}: {'ALIVE' if alive else 'DEAD'}")
    if not alive:
        problems.append(f"pin in {p} is dead; this record asserts pins_touched=0")

rec = json.loads(Path(
    "config/artifacts/incidents/mpralegnet_deliberate_non_changes_20260825.json"
).read_text(encoding="utf-8"))
assert rec["files_modified"] == 0, "record claims files_modified=0"
assert rec["pins_touched"] == 0, "record claims pins_touched=0"
assert {e["site"] for e in rec["non_changes"]} == {"D", "E", "F"}
print(f"\nrecord declares {len(rec['non_changes'])} deliberate non-changes: D, E, F")
if problems:
    raise SystemExit("record is stale:\n" + "\n".join(problems))
print("deliberate-non-changes record re-derives; both pins alive; base untouched")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_mpralegnet_deliberate_non_changes \
  tests.unit.test_mpralegnet_terminal_disposition \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
tail -5 "${STAGING}/validation/unittest.txt"
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_registry_files | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp config/artifacts/incidents/mpralegnet_deliberate_non_changes_20260825.json "${STAGING}/audit/"
sha256sum "${ROOT}/config/artifacts/incidents/mpralegnet_deliberate_non_changes_20260825.json" \
          "${ROOT}/tests/unit/test_mpralegnet_deliberate_non_changes.py" \
          "${ROOT}/slurm/freeze_mpralegnet_deliberate_non_changes.sbatch" > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"mpralegnet_deliberate_non_changes_record",
    "files_modified":0,"pins_touched":0,"cascade_started":False,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
