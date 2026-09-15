#!/usr/bin/env bash
#SBATCH --job-name=model-check-920
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
TARGET="${ROOT}/executions/pin-companion-refresh-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_refresh_pin_companion_metadata \
  > "${STAGING}/validation/unittest_companion.txt" 2>&1

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/refresh_pin_companion_metadata.py \
  --root "${ROOT}" --output "${STAGING}/audit/companion" --apply \
  --record config/artifacts/models/evidence_count/checkpoints.json \
  --record config/artifacts/models/evidence_graph_native_baselines/checkpoints.json \
  --record config/artifacts/models/late_fusion/checkpoints.json \
  --record config/artifacts/models/nearest_gene/checkpoints.json \
  --record config/artifacts/models/negative_controls/checkpoints.json \
  --record config/artifacts/models/cell_bulk_classical_baselines/checkpoints.json \
  --record config/artifacts/models/cell_bulk_classical_baselines/development_crosswalk.json \
  > "${STAGING}/validation/companion.stdout.txt" 2>&1

PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/post_state.txt"
import json
from pathlib import Path
from scripts.refresh_pin_companion_metadata import stale_companions
root=Path(".")
recs=json.loads(Path([p for p in Path(".").glob("executions/pin-companion-refresh-*.staging/audit/companion/receipt.json")][0]).read_text())
print(f"stale companions found: {recs['total_stale']}")
left=0
for r in recs["records"]:
    n=stale_companions(root, root/r)
    left+=len(n)
    for c in recs["records"][r]["changes"]:
        print(f"   {r.split('/')[-2]}/{r.split('/')[-1]}{c['json_path']}.{c['key']}  {c['from']} -> {c['to']}")
print(f"stale companions remaining in these records: {left}")
assert left==0, "companions still stale"
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_refresh_pin_companion_metadata tests.unit.test_refresh_authorized_pins \
  tests.unit.test_check_contract_suite_baseline tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_refresh_pin_companion_metadata --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
sha256sum "${ROOT}/scripts/refresh_pin_companion_metadata.py" > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"pin_companion_metadata_refresh","records":7,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
