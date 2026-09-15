#!/usr/bin/env bash
#SBATCH --job-name=model-check-993
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
CONTRACT="${ROOT}/config/campaigns/scbasset_five_seed_rectangle_20260825.json"
TARGET="${ROOT}/executions/scbasset-five-seed-rectangle-admission-20260825"
STAGE="${TARGET}.staging-${SLURM_JOB_ID}"
FAILED="${TARGET}.failed-${SLURM_JOB_ID}"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite scBasset five-seed admission: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "3b7d6ef7c985d126fd0fa16bb7b5c51594e38309cbeb0c4774f03444b5dc809b config/campaigns/scbasset_five_seed_rectangle_20260825.json" \
  "5b0cdcc156173e801be557df83e03a35952cb283694790b14c2d44101a78dac2 config/artifacts/incidents/scbasset_broad_rg_prediction_line_exposure_20260825.json" \
  "69970083193e7566ce8526a3defd1831ac169556be78998c6fd000c7c897a1bc scripts/audit_scbasset_five_seed_rectangle_admission.py" \
  "7db96ab795a7e83373fe19fc0b9dbbf55528c1acb5d0f1cf4761cdb9386a903b tests/unit/test_audit_scbasset_five_seed_rectangle_admission.py" \
  "810be9dc07bef8f8c8a6d011641482848628f32cc002754dfeb1cc160676769b scripts/validate_scbasset_broad_rg_prediction_line_exposure.py" \
  "0fa86b9564e2452fd5f450118e6ba1adbf3731ee910796cb2a545eaa957c818d tests/unit/test_validate_scbasset_broad_rg_prediction_line_exposure.py" \
  "08402a3900884386417951bb84e8267addd3a8d8181984f2532af59e6fe67bc1 slurm/run_scbasset_five_seed_extension_bundle.sbatch" \
  "77d98afcbfc46d176f81a5c2a1268e3237d8a4b01f612c1e3060ee8fc5883097 config/campaigns/gpu_bundle_queue/064-model-training-605.json" \
  "68f3eb8e8c6e2cfed15f03b4b43ac0eadb444c097720e1dc7ad55c4c881ee75e slurm/scbasset_train_donor0_genomic0_seed11.sbatch" \
  "c647db4f41a0f1ee1260de68c6fc476f663a21d65be95b59c19cc062578cd8d5 slurm/scbasset_predict_valid_donor0_genomic0_seed11.sbatch"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_scbasset_five_seed_rectangle_admission \
  tests.unit.test_validate_scbasset_broad_rg_prediction_line_exposure \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/validate_scbasset_broad_rg_prediction_line_exposure.py \
  --incident "${ROOT}/config/artifacts/incidents/scbasset_broad_rg_prediction_line_exposure_20260825.json" \
  --output "${STAGE}/incident/validation_receipt.json" \
  > "${STAGE}/validation/incident.stdout.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_scbasset_five_seed_rectangle_admission.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/admission/receipt.json" \
  > "${STAGE}/validation/audit.stdout.txt"

SCBASSET_FIVE_SEED_AUTHORITY=scbasset_five_seed_rectangle_v1 \
SCBASSET_FIVE_SEED_CONTRACT="${CONTRACT}" \
SCBASSET_FIVE_SEED_CONTRACT_SHA256=3b7d6ef7c985d126fd0fa16bb7b5c51594e38309cbeb0c4774f03444b5dc809b \
SCBASSET_PREFLIGHT_ONLY=1 \
SCBASSET_PREFLIGHT_STAGE="${STAGE}/validation/bundle-preflight" \
bash slurm/run_scbasset_five_seed_extension_bundle.sbatch \
  > "${STAGE}/validation/bundle-preflight.stdout.txt" \
  2> "${STAGE}/validation/bundle-preflight.stderr.txt"

python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/campaigns/scbasset_five_seed_rectangle_20260825.json \
  config/artifacts/incidents/scbasset_broad_rg_prediction_line_exposure_20260825.json \
  config/campaigns/gpu_bundle_queue/064-model-training-605.json \
  scripts/audit_scbasset_five_seed_rectangle_admission.py \
  tests/unit/test_audit_scbasset_five_seed_rectangle_admission.py \
  scripts/validate_scbasset_broad_rg_prediction_line_exposure.py \
  tests/unit/test_validate_scbasset_broad_rg_prediction_line_exposure.py \
  slurm/run_scbasset_five_seed_extension_bundle.sbatch \
  slurm/scbasset_train_donor0_genomic0_seed11.sbatch \
  slurm/scbasset_predict_valid_donor0_genomic0_seed11.sbatch \
  slurm/admit_scbasset_five_seed_rectangle_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${SLURM_JOB_ID}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
job_id = sys.argv[2]
receipt = json.loads((root / "admission/receipt.json").read_text(encoding="utf-8"))
incident = json.loads(
    (root / "incident/validation_receipt.json").read_text(encoding="utf-8")
)
preflight = json.loads(
    (root / "validation/bundle-preflight/wrapper_preflight.json").read_text(
        encoding="utf-8"
    )
)
fits = list(
    csv.DictReader(
        (root / "validation/bundle-preflight/fit_plan.tsv").open(), delimiter="\t"
    )
)
stderr = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if (
    receipt.get("status")
    != "pass_exact_21_missing_fit_rectangle_covered_by_two_central_bundles"
    or receipt.get("expected_fits") != 25
    or receipt.get("complete_fits_at_freeze") != 4
    or receipt.get("missing_fits_admitted") != 21
    or receipt.get("expected_valid_prediction_views") != 25
    or receipt.get("prediction_view")
    != "valid_training_cell_state_mean_sequence_profile"
    or receipt.get("donor_context") != "none_sequence_only"
    or receipt.get("observed_query_atac_consumed") is not False
    or receipt.get("rna_context_consumed") is not False
    or receipt.get("observed_multiome_inductive_claim_allowed") is not False
    or receipt.get("existing_bundle_id") != "model-training-604"
    or receipt.get("existing_bundle_missing_fits") != 11
    or receipt.get("extension_bundle_id") != "model-training-605"
    or receipt.get("extension_bundle_missing_fits") != 10
    or receipt.get("central_bundle_count") != 2
    or receipt.get("logical_tasks") != 42
    or receipt.get("requested_new_gpu_allocations") != 1
    or receipt.get("central_queue_registration_present") is not True
    or receipt.get("central_queue_unclaimed_during_validation") is not True
    or receipt.get("direct_gpu_submission_authorized") is not False
    or receipt.get("partial_ranking_authorized") is not False
    or receipt.get("development_evaluator_open") is not False
    or receipt.get("admission_audit_scbasset_prediction_values_read") is not False
    or receipt.get("metric_values_read") is not False
    or receipt.get("raw_outcomes_read") is not False
    or receipt.get("sealed_labels_or_outcomes_read") is not False
    or receipt.get("development_search_incident_recorded") is not True
    or receipt.get("development_search_incident_id")
    != "scbasset_broad_rg_prediction_line_exposure_20260825"
    or receipt.get("development_search_incident_sha256")
    != "5b0cdcc156173e801be557df83e03a35952cb283694790b14c2d44101a78dac2"
    or receipt.get("development_search_incident_structure_validated") is not True
    or receipt.get("unrelated_prediction_lines_accidentally_emitted") is not True
    or receipt.get("incident_values_used") is not False
    or receipt.get("gpu_job_submitted_by_admission") is not False
    or preflight.get("status") != "pass"
    or preflight.get("bundle_id") != "model-training-605"
    or preflight.get("fits") != 10
    or preflight.get("logical_tasks") != 20
    or preflight.get("seed_gate_change_only") is not True
    or preflight.get("direct_gpu_submission_performed") is not False
    or incident.get("status") != "pass_structure_only_incident_validation"
    or incident.get("incident_id")
    != "scbasset_broad_rg_prediction_line_exposure_20260825"
    or incident.get("incident_sha256")
    != "5b0cdcc156173e801be557df83e03a35952cb283694790b14c2d44101a78dac2"
    or incident.get("incident_values_used") is not False
    or incident.get("referenced_prediction_files_reopened_by_validator") is not False
    or len(fits) != 10
    or {(int(row["outer_fold"]), int(row["seed"])) for row in fits}
    != {(fold, seed) for fold in range(5) for seed in (20260827, 20260828)}
    or "Ran 9 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("scBasset five-seed admission receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "scbasset_five_seed_rectangle_admission",
        "model_id": "scbasset",
        "dataset_id": "gse296875",
        "slurm_job_id": job_id,
        "tests": 9,
        "expected_fits": 25,
        "complete_fits_at_freeze": 4,
        "missing_fits_admitted": 21,
        "expected_valid_prediction_views": 25,
        "prediction_view": "valid_training_cell_state_mean_sequence_profile",
        "donor_context": "none_sequence_only",
        "observed_query_atac_consumed": False,
        "rna_context_consumed": False,
        "observed_multiome_inductive_claim_allowed": False,
        "existing_bundle_id": "model-training-604",
        "existing_bundle_missing_fits": 11,
        "extension_bundle_id": "model-training-605",
        "extension_bundle_missing_fits": 10,
        "central_bundle_count": 2,
        "logical_tasks": 42,
        "requested_new_gpu_allocations": 1,
        "requested_new_gpu_hours": 12,
        "queue_gate_materialized_by_this_artifact": True,
        "direct_gpu_submission_performed": False,
        "partial_ranking_performed": False,
        "admission_audit_scbasset_prediction_values_read": False,
        "metric_values_read": False,
        "raw_outcomes_read": False,
        "sealed_labels_or_outcomes_read": False,
        "development_search_incident_recorded": True,
        "development_search_incident_id": incident["incident_id"],
        "development_search_incident_sha256": incident["incident_sha256"],
        "development_search_incident_structure_validated": True,
        "unrelated_prediction_lines_accidentally_emitted": True,
        "incident_values_used": False,
        "contract_sha256": receipt["contract_sha256"],
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'scBasset five-seed rectangle admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
