#!/usr/bin/env bash
#SBATCH --job-name=model-validate-206
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
TARGET="${ROOT}/executions/model-validate-206-${SLURM_JOB_ID}"
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

python3 -m unittest -v tests.unit.test_audit_bulk_protein_transportability \
  > "${STAGE}/validation/unit-test.stdout.txt" \
  2> "${STAGE}/validation/unit-test.stderr.txt"
python3 scripts/audit_bulk_protein_transportability.py \
  --project-root "${ROOT}" \
  --output "${STAGE}/audit"
python3 -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum scripts/audit_bulk_protein_transportability.py \
  tests/unit/test_audit_bulk_protein_transportability.py \
  config/artifacts/model-authorities/bulk_protein_transportability/census.json \
  slurm/validate_model_lane_206_cpu.sbatch \
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
    (stage / "audit/bulk_protein_transportability_audit.json").read_text(
        encoding="utf-8"
    )
)
if (
    audit["status"] != "pass_fail_closed"
    or audit["cohort_count"] != 10
    or audit["next_outcome_blind_campaign"] != "gse268273_raw_reprocessing_v1"
    or audit["additional_parallel_campaign_justified"]
    or audit["submission_authorized"]
    or audit["models_fit"]
    or audit["models_scored"]
    or audit["sealed_sources_accessed"]
    or audit["bulk_champion_created"]
    or audit["protein_champion_created"]
    or audit["resource_authorities_modified"]
):
    raise SystemExit("bulk/protein transportability audit disposition differs")
freeze_tree(
    stage,
    metadata={
        "artifact_class": "bulk_protein_transportability_reconciliation",
        "cohort_count": 10,
        "next_outcome_blind_campaign": "gse268273_raw_reprocessing_v1",
        "new_jobs_created": 0,
        "sealed_sources_accessed": [],
        "status": "passed_fail_closed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
