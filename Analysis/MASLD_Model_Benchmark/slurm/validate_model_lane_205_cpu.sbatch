#!/usr/bin/env bash
#SBATCH --job-name=model-validate-205
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
TARGET="${ROOT}/executions/model-validate-205-${SLURM_JOB_ID}"
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

python3 -m unittest -v tests.unit.test_audit_remaining_cell_specialists \
  > "${STAGE}/validation/unit-test.stdout.txt" \
  2> "${STAGE}/validation/unit-test.stderr.txt"
python3 scripts/audit_remaining_cell_specialists.py \
  --project-root "${ROOT}" \
  --output "${STAGE}/audit"
python3 -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum scripts/audit_remaining_cell_specialists.py \
  tests/unit/test_audit_remaining_cell_specialists.py \
  config/artifacts/model-authorities/remaining_cell_specialists/census.json \
  slurm/validate_model_lane_205_cpu.sbatch \
  > "${STAGE}/source.sha256"

python3 - "${ROOT}" "${STAGE}" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
stage = Path(sys.argv[2])
sys.path.insert(0, str(root / "src"))
from masld_bench.artifacts import freeze_tree

audit = json.loads(
    (stage / "audit/cell_specialist_audit.json").read_text(encoding="utf-8")
)
if (
    audit["status"] != "pass_fail_closed"
    or audit["next_executable_50k_model"] != "scimilarity_v1_1"
    or audit["submission_authorized"]
    or audit["models_trained"]
    or audit["models_scored"]
    or audit["sealed_sources_accessed"]
    or audit["sealed_champion_created"]
):
    raise SystemExit("remaining cell-specialist audit disposition differs")
freeze_tree(
    stage,
    metadata={
        "artifact_class": "remaining_cell_specialist_reconciliation",
        "models": audit["models_reconciled"],
        "next_executable_50k_model": "scimilarity_v1_1",
        "new_jobs_created": 0,
        "sealed_sources_accessed": [],
        "status": "passed_fail_closed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
