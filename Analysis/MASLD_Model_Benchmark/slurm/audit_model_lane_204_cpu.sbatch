#!/usr/bin/env bash
#SBATCH --job-name=model-audit-204
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
TARGET="${ROOT}/executions/model-audit-204-${SLURM_JOB_ID}"
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

python3 -m unittest -v tests.unit.test_audit_perturbation_specialist_census \
  > "${STAGE}/validation/unit-test.stdout.txt" \
  2> "${STAGE}/validation/unit-test.stderr.txt"
python3 scripts/audit_perturbation_specialist_census.py \
  --project-root "${ROOT}" \
  --output "${STAGE}/audit"
python3 -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum scripts/audit_perturbation_specialist_census.py \
  tests/unit/test_audit_perturbation_specialist_census.py \
  config/artifacts/model-authorities/perturbation_specialists/census.json \
  slurm/audit_model_lane_204_cpu.sbatch \
  > "${STAGE}/source.sha256"

python3 - "${ROOT}" "${STAGE}" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
stage = Path(sys.argv[2])
sys.path.insert(0, str(root / "src"))
from masld_bench.artifacts import freeze_tree

audit = json.loads((stage / "audit/census_audit.json").read_text(encoding="utf-8"))
if (
    audit["status"] != "pass_fail_closed"
    or audit["sealed_sources"]
    or audit["champion_eligible"]
    or audit["production_jobs_authorized"]
    or not audit["native_subendpoints_preserved"]
    or not audit["biological_replication_preserved"]
):
    raise SystemExit("perturbation specialist census disposition differs")
freeze_tree(
    stage,
    metadata={
        "artifact_class": "perturbation_specialist_census_audit",
        "models": audit["models_reconciled"],
        "new_jobs_created": 0,
        "sealed_sources": [],
        "champion_eligible": False,
        "status": "passed_fail_closed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
