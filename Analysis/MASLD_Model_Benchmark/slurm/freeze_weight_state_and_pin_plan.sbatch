#!/usr/bin/env bash
#SBATCH --job-name=model-check-906
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
TARGET="${ROOT}/executions/weight-state-and-pin-plan-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# Re-derive the scgpt and uce weight digests from the bodies themselves.
# The acquisition SHA256SUMS receipt is corroboration, not the authority.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/weight_digest_rederivation.txt"
import hashlib, json
from pathlib import Path
CASES = [
    ("scgpt", "executions/scgpt-acquisition-21015743/scgpt_continual/best_model.pt",
     207861754, "ad0252a1971e0cd619b7116dbab3177432236c4537225d54280a2aa7e5fe402a"),
    ("scgpt", "executions/scgpt-acquisition-21015743/scgpt_whole_human/best_model.pt",
     205385258, "6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922"),
    ("uce", "executions/uce-acquisition-21015745/33l_8ep_1024t_1280.torch",
     5686283829, "3f458726196308e171611ed28394b55865708749f82af68cfd7771d5dfee661e"),
    ("uce", "executions/uce-acquisition-21015745/4layer_model.torch",
     3403514339, "acb28f3f0a1d803e4a4ffe891b9bab38bf93c84762dc06b2452f0d515da91560"),
]
bad = []
for model, rel, size, expected in CASES:
    p = Path(rel)
    if not p.is_file():
        bad.append(f"{rel} missing"); continue
    actual_size = p.stat().st_size
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(64 * 1024 * 1024), b""):
            h.update(blk)
    d = h.hexdigest()
    ok = (d == expected and actual_size == size)
    print(f"{'ok ' if ok else 'BAD'} {model:6s} {actual_size:>13,}  {d}  {rel}")
    if not ok: bad.append(rel)
if bad: raise SystemExit("weight digest re-derivation failed: " + ", ".join(bad))

sc = json.loads(Path("config/artifacts/models/scgpt/checkpoints.json").read_text())
uc = json.loads(Path("config/artifacts/models/uce/checkpoints.json").read_text())
print(f"\nscgpt bundle weight_content_downloaded = {sc['weight_content_downloaded']} "
      f"(bodies present -> should be True), weight_license = {sc['weight_license']}")
print(f"uce   bundle weight_content_downloaded = {uc['weight_content_downloaded']} "
      f"(bodies present -> correct), weight_license = {uc['weight_license']}")
assert sc["weight_content_downloaded"] is False, "scgpt bundle changed; re-audit"
assert uc["weight_content_downloaded"] is True, "uce bundle changed; re-audit"
print("investigation findings still hold; no file was modified")
PY

# Re-derive the pin plan against the live tree so a stale plan cannot be frozen.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/pin_plan_rederivation.txt"
import hashlib, json
from pathlib import Path
plan = json.loads(Path("config/campaigns/stale_pin_refresh_plan_20260825.json").read_text())
assert plan["execution_authorized"] is False and plan["status"] == "PLAN_ONLY_NOT_EXECUTED"
bad = []
for s in plan["moved_sources"]:
    p = Path(s["path"])
    live = hashlib.sha256(p.read_bytes()).hexdigest()
    if live != s["live_digest"]:
        bad.append(f"{s['path']}: plan says {s['live_digest'][:12]} but is {live[:12]}")
print(f"moved sources re-derived: {len(plan['moved_sources'])}, mismatches: {len(bad)}")
for b in bad: print("  ", b)
if bad: raise SystemExit("plan is stale against the live tree")
print(f"dead pins: {len(plan['dead_pins'])}, cascade waves: {len(plan['cascade_waves'])}")
print("pin refresh plan re-derives")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_campaign_census_scoping \
  tests.unit.test_check_contract_suite_baseline \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp config/artifacts/incidents/scgpt_uce_weight_state_investigation_20260825.json "${STAGING}/audit/"
cp config/campaigns/stale_pin_refresh_plan_20260825.json "${STAGING}/audit/"
sha256sum \
  "${ROOT}/config/artifacts/incidents/scgpt_uce_weight_state_investigation_20260825.json" \
  "${ROOT}/config/campaigns/stale_pin_refresh_plan_20260825.json" \
  "${ROOT}/config/artifacts/models/scgpt/checkpoints.json" \
  "${ROOT}/config/artifacts/models/uce/checkpoints.json" \
  "${ROOT}/slurm/freeze_weight_state_and_pin_plan.sbatch" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class": "weight_state_investigation_and_pin_refresh_plan",
    "investigation_scope": "report_only", "plan_scope": "plan_only_not_executed",
    "files_modified": 0, "checkpoint_opened": False, "weights_deserialized": False,
    "sealed_outcomes_read": False, "models_fit": False, "models_scored": False,
}))
PY
mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
