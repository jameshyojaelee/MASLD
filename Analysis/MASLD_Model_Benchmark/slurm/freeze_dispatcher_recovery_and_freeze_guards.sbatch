#!/usr/bin/env bash
#SBATCH --job-name=model-check-925
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
TARGET="${ROOT}/executions/dispatcher-recovery-and-freeze-guards-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"; FAILED="${TARGET}.failed"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
PROMOTION_PENDING=0
preserve_failed() {
  if [[ -d "${STAGING}" && ! -e "${FAILED}" ]]; then mv "${STAGING}" "${FAILED}"; fi
  if [[ "${PROMOTION_PENDING}" -eq 1 && -e "${TARGET}" && ! -e "${FAILED}" ]]; then mv "${TARGET}" "${FAILED}"; fi
}
trap preserve_failed EXIT
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# The guard's intent is "am I changing the LIVE dispatcher", not "is any
# dispatcher busy". Those differ: the second blocks unrelated work whenever the
# dispatcher is running, and would permit a live edit whenever it happened to be
# idle. Test the property the guard is actually about.
#
# The live entry point is scripts/gpu_dependency_dispatcher.py. Job 21100574
# runs it and it must not change underneath that process. This campaign fixed
# only the LEGACY scripts/gpu_bundle_dispatcher.py, whose dispatch_once the live
# path never calls, so freezing that is safe while the dispatcher runs.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PYGUARD' | tee "${STAGING}/audit/live_dispatcher_untouched.txt"
import hashlib
from pathlib import Path
LIVE = Path("scripts/gpu_dependency_dispatcher.py")
LEGACY = Path("scripts/gpu_bundle_dispatcher.py")
print(f"live entry point {LIVE} sha256 {hashlib.sha256(LIVE.read_bytes()).hexdigest()[:16]}")
print(f"legacy module    {LEGACY} (fixed, inert)")
text = LIVE.read_text(encoding="utf-8")
assert '.exists() or summary.get("enabled") is False' in text, \
    "the live dispatcher changed; re-audit before freezing anything about it"
assert '"claimed_before_submission"' in text
print("live dispatcher still carries the unfixed defect at :524 / :637, as recorded")
assert "legacy.dispatch_once" not in text, \
    "the live path calls legacy.dispatch_once; the legacy fix is NOT inert"
print("live path does not call legacy.dispatch_once; legacy fix confirmed inert")
PYGUARD

# Freeze-guard ordering, asserted rather than eyeballed.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/freeze_guard_ordering.txt"
import re
from pathlib import Path
bad=[]
for name in ("build_gse256398_qc_sensitivity","build_gse256398_transcriptformer_shards"):
    text=Path(f"slurm/{name}.sbatch").read_text().splitlines()
    def line_of(pat):
        return next((i for i,l in enumerate(text,1) if re.search(pat,l)), None)
    trap  = line_of(r"^trap preserve_failed")
    mv    = line_of(r'^mv "\$\{STAGE\}" "\$\{TARGET\}"')
    setp  = line_of(r"^PROMOTION_PENDING=1")
    clear = line_of(r"^PROMOTION_PENDING=0$")
    untrap= line_of(r"^trap - EXIT")
    clears=[i for i,l in enumerate(text,1) if l.strip()=="PROMOTION_PENDING=0"]
    final_clear=max(clears)
    print(f"{name}: trap@{trap} mv@{mv} pending@{setp} clear@{final_clear} untrap@{untrap}")
    if not (trap < mv < setp < final_clear < untrap):
        bad.append(f"{name}: ordering wrong")
    # the trap must actually cover TARGET, not only STAGE
    body="\n".join(text)
    if "PROMOTION_PENDING}\" -eq 1" not in body or 'mv "${TARGET}" "${FAILED}"' not in body:
        bad.append(f"{name}: trap does not rename TARGET")
if bad: raise SystemExit("freeze-guard ordering wrong:\n"+"\n".join(bad))
print("both freeze guards: trap installed before mv, pending set after mv, cleared only after the checks")
PY

set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_gpu_dispatcher_claim_recovery tests.unit.test_config_write_transaction \
  tests.unit.test_check_contract_suite_baseline tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_gpu_dispatcher_claim_recovery --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
cp config/artifacts/incidents/verification_method_register_20260825.json "${STAGING}/audit/"
sha256sum "${ROOT}/scripts/gpu_bundle_dispatcher.py" \
  "${ROOT}/slurm/build_gse256398_qc_sensitivity.sbatch" \
  "${ROOT}/slurm/build_gse256398_transcriptformer_shards.sbatch" \
  "${ROOT}/config/artifacts/incidents/verification_method_register_20260825.json" \
  > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"dispatcher_recovery_and_freeze_guard_fixes",
    "findings_fixed":[2,4],"finding_1_dropped":True,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"
PROMOTION_PENDING=1
sha256sum "${TARGET}/ARTIFACTS.json"
PROMOTION_PENDING=0
trap - EXIT
echo "FROZEN: ${TARGET}"
