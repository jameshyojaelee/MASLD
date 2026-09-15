#!/usr/bin/env bash
#SBATCH --job-name=model-check-982
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT="${ROOT}/config/artifacts/models/chrombpnet/full_head_production_readiness_20260825.json"
TARGET="${ROOT}/executions/model-check-982-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite ChromBPNet readiness validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "cb7c2173567b9ba1e1018615d932225a7c0362ac2310c7025daf5df2254d813d config/artifacts/models/chrombpnet/full_head_production_readiness_20260825.json" \
  "4df6fa403b500cc416d1f336995c711a65a0b8b5e673f21701feb213ae342a45 scripts/audit_chrombpnet_full_head_readiness.py" \
  "0f621aec0805e0282b9835f8fe0a2570e2f753ef6fb1ba7d6e13132093569ea4 tests/unit/test_audit_chrombpnet_full_head_readiness.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_chrombpnet_full_head_readiness \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_chrombpnet_full_head_readiness.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/readiness/receipt.json" \
  > "${STAGE}/validation/audit.stdout.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/artifacts/models/chrombpnet/full_head_production_readiness_20260825.json \
  scripts/audit_chrombpnet_full_head_readiness.py \
  tests/unit/test_audit_chrombpnet_full_head_readiness.py \
  slurm/validate_chrombpnet_full_head_readiness_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "readiness/receipt.json").read_text(encoding="utf-8"))
stderr = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if (
    receipt.get("status")
    != "pass_full_head_reconciled_waiting_for_rectangle_and_license"
    or receipt.get("primary_biological_view") != "chrombpnet_nobias"
    or receipt.get("assay_qc_view") != "chrombpnet_full"
    or receipt.get("chrombpnet_complete_fits_at_freeze") != 19
    or receipt.get("chrombpnet_missing_fits_at_freeze") != 106
    or receipt.get("full_compatible_rectangle") is not False
    or receipt.get("license_authority_contradiction_present") is not True
    or receipt.get("weight_redistribution_authorized") is not False
    or receipt.get("open_champion_eligible_now") is not False
    or receipt.get("valid_next_gpu_unit_already_registered") is not True
    or receipt.get("next_dispatcher_bundle_id") != "model-training-801"
    or receipt.get("direct_gpu_submission_authorized") is not False
    or receipt.get("prediction_values_read") is not False
    or receipt.get("development_outcomes_read") is not False
    or receipt.get("sealed_labels_or_outcomes_read") is not False
    or receipt.get("gpu_job_submitted") is not False
    or "Ran 5 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("ChromBPNet full-head readiness receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "chrombpnet_full_head_production_readiness",
        "model_id": "chrombpnet",
        "artifact_role": "outcome_blind_head_rectangle_license_and_dispatch_reconciliation",
        "tests": 5,
        "primary_biological_view": "chrombpnet_nobias",
        "assay_qc_view": "chrombpnet_full",
        "complete_fits_at_freeze": 19,
        "missing_fits_at_freeze": 106,
        "full_compatible_rectangle": False,
        "license_authority_contradiction_present": True,
        "weight_redistribution_authorized": False,
        "open_champion_eligible_now": False,
        "valid_next_gpu_unit_already_registered": True,
        "next_dispatcher_bundle_id": "model-training-801",
        "gpu_job_submitted": False,
        "prediction_values_read": False,
        "development_outcomes_read": False,
        "sealed_labels_or_outcomes_read": False,
        "contract_sha256": receipt["contract_sha256"],
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet full-head readiness validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
