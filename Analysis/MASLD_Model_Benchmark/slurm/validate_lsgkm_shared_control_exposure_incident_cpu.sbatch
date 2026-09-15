#!/usr/bin/env bash
#SBATCH --job-name=model-check-619
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
INCIDENT="${ROOT}/config/artifacts/incidents/lsgkm_shared_control_prediction_exposure_20260825.json"
TARGET="${ROOT}/executions/model-check-619-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite incident validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' 8485d619d45f9a9d22297c179193364619a205f90f49efc0d4e6db26c494c314 "${INCIDENT}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INCIDENT}" "${STAGE}" <<'PY'
from pathlib import Path
import hashlib
import json
import sys
from masld_bench.artifacts import freeze_tree

incident_path = Path(sys.argv[1]).resolve(strict=True)
output = Path(sys.argv[2]).resolve(strict=True)
incident = json.loads(incident_path.read_text(encoding="utf-8"))
expected_false = (
    "reporter_outcomes_opened",
    "reporter_counts_opened",
    "development_labels_opened",
    "sealed_features_opened",
    "sealed_labels_or_outcomes_opened",
    "benchmark_metrics_opened_or_calculated",
    "additional_control_prediction_values_opened",
    "prediction_value_used_for_model_training",
    "prediction_value_used_for_model_selection",
    "prediction_value_used_for_threshold_or_calibration_choice",
)
if (
    incident.get("schema_version") != "masld-bench-protocol-exposure-incident-v1"
    or incident.get("incident_id") != "lsgkm_shared_control_prediction_exposure_20260825"
    or incident.get("exposure_scope", {}).get("data_rows_opened") != 1
    or incident.get("exposure_scope", {}).get("prediction_values_opened") != 1
    or incident.get("exposure_scope", {}).get("exposed_head_id") != "zero"
    or incident.get("impact_assessment", {}).get("control_value_blinding_breach") is not True
    or any(incident.get("explicit_exclusions", {}).get(field) is not False for field in expected_false)
    or incident.get("containment", {}).get("incident_artifact_contains_no_prediction_value") is not True
    or not all(incident.get("required_downstream_handling", {}).values())
):
    raise SystemExit("incident structure or containment contract differs")
digest = hashlib.sha256(incident_path.read_bytes()).hexdigest()
receipt = {
    "schema_version": "masld-bench-protocol-exposure-incident-validation-v1",
    "status": "pass_structure_only_incident_validation",
    "incident_id": incident["incident_id"],
    "incident_sha256": digest,
    "referenced_prediction_member_opened": False,
    "additional_prediction_values_opened": False,
    "outcomes_opened": False,
    "metrics_calculated": False,
}
(output / "validation_receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
freeze_tree(output, {
    "artifact_class": "protocol_exposure_incident_validation",
    **receipt,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Incident validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
