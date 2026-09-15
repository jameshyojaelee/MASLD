#!/usr/bin/env bash
#SBATCH --job-name=model-audit-202
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
module load Python/3.11.3-GCCcore-12.3.0

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/model-audit-202-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
test ! -e "${TARGET}"
test ! -e "${STAGE}"
test ! -e "${FAILED}"
mkdir -p "${STAGE}/validation"
export PYTHONPATH="${ROOT}/src"
cd "${ROOT}"

python3 -m unittest -v tests.unit.test_audit_enformer_sei_native_lanes \
  > "${STAGE}/validation/unit-test.stdout.txt" \
  2> "${STAGE}/validation/unit-test.stderr.txt"
python3 scripts/audit_enformer_sei_native_lanes.py \
  --project-root "${ROOT}" \
  --output "${STAGE}/audit"
python3 -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum scripts/audit_enformer_sei_native_lanes.py \
  tests/unit/test_audit_enformer_sei_native_lanes.py \
  slurm/audit_model_lane_202_cpu.sbatch \
  > "${STAGE}/source.sha256"

python3 - "${ROOT}" "${STAGE}" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
stage = Path(sys.argv[2])
sys.path.insert(0, str(root / "src"))
from masld_bench.artifacts import freeze_tree

audit = json.loads((stage / "audit/lane_audit.json").read_text(encoding="utf-8"))
if (
    audit["status"] != "pass"
    or audit["enformer"]["native_sonnet_parity_established"]
    or audit["enformer"]["open_champion_eligible"]
    or audit["sei"]["open_champion_eligible"]
    or audit["queue_disposition"]["new_gpu_job_required"]
    or audit["project_data_read_by_audit"]
    or audit["outcomes_read_by_audit"]
    or audit["sealed_labels_read_by_audit"]
):
    raise SystemExit("Enformer/Sei audit disposition differs")
freeze_tree(
    stage,
    metadata={
        "artifact_class": "enformer_sei_native_lane_audit",
        "models": ["enformer_crested_restricted_port", "sei"],
        "new_gpu_jobs_created": 0,
        "outcomes_read": False,
        "sealed_labels_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
