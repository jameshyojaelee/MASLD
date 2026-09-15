#!/usr/bin/env bash
#SBATCH --job-name=model-check-988
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
CONTRACT="${ROOT}/config/campaigns/scbasset_three_seed_rectangle_20260825.json"
TARGET="${ROOT}/executions/scbasset-three-seed-rectangle-admission-20260825"
STAGE="${TARGET}.staging-${SLURM_JOB_ID}"
FAILED="${TARGET}.failed-${SLURM_JOB_ID}"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite scBasset rectangle admission: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "66b62badd3e8afa61815afdfaa18ce81100f39e6dfd0918094c8127e3c09a695 config/campaigns/scbasset_three_seed_rectangle_20260825.json" \
  "007d623f8b6fda28c81de1790006a17cee8db0831fce35454967dc8381eeff4e scripts/audit_scbasset_three_seed_rectangle_admission.py" \
  "5338ee9cae2e458a1b9d85f85b3d63a0ee8e71d74761e9323a75423568c6cc7a tests/unit/test_audit_scbasset_three_seed_rectangle_admission.py" \
  "ac9d3e459c17b998e0f314797f95a705e79bd69596022b7d3d7f03fd0639d549 slurm/run_scbasset_three_seed_rectangle_bundle.sbatch" \
  "68f3eb8e8c6e2cfed15f03b4b43ac0eadb444c097720e1dc7ad55c4c881ee75e slurm/scbasset_train_donor0_genomic0_seed11.sbatch" \
  "c647db4f41a0f1ee1260de68c6fc476f663a21d65be95b59c19cc062578cd8d5 slurm/scbasset_predict_valid_donor0_genomic0_seed11.sbatch" \
  "f2f66484704ac5e3719a5101980f018e6492d3eebc2a8acba99240a07ac4ff3c config/campaigns/gpu_bundle_queue/063-model-training-604.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest -v \
  tests.unit.test_audit_scbasset_three_seed_rectangle_admission \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_scbasset_three_seed_rectangle_admission.py \
  --project-root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}/admission/receipt.json" \
  > "${STAGE}/validation/audit.stdout.txt"

SCBASSET_RECTANGLE_AUTHORITY=scbasset_three_seed_rectangle_v1 \
SCBASSET_RECTANGLE_CONTRACT="${CONTRACT}" \
SCBASSET_RECTANGLE_CONTRACT_SHA256=66b62badd3e8afa61815afdfaa18ce81100f39e6dfd0918094c8127e3c09a695 \
SCBASSET_PREFLIGHT_ONLY=1 \
SCBASSET_PREFLIGHT_STAGE="${STAGE}/validation/bundle-preflight" \
bash slurm/run_scbasset_three_seed_rectangle_bundle.sbatch \
  > "${STAGE}/validation/bundle-preflight.stdout.txt" \
  2> "${STAGE}/validation/bundle-preflight.stderr.txt"

python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/campaigns/scbasset_three_seed_rectangle_20260825.json \
  config/campaigns/gpu_bundle_queue/063-model-training-604.json \
  scripts/audit_scbasset_three_seed_rectangle_admission.py \
  tests/unit/test_audit_scbasset_three_seed_rectangle_admission.py \
  slurm/run_scbasset_three_seed_rectangle_bundle.sbatch \
  slurm/scbasset_train_donor0_genomic0_seed11.sbatch \
  slurm/scbasset_predict_valid_donor0_genomic0_seed11.sbatch \
  slurm/admit_scbasset_three_seed_rectangle_cpu.sbatch \
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
    != "pass_exact_11_fit_rectangle_admitted_for_central_dispatch"
    or receipt.get("complete_fits_at_freeze") != 4
    or receipt.get("missing_fits_admitted") != 11
    or receipt.get("fixed_seeds") != [20260824, 20260825, 20260826]
    or receipt.get("legacy_seeds_excluded") != [11, 29, 47, 71, 101]
    or receipt.get("central_bundle_id") != "model-training-604"
    or receipt.get("central_bundle_count") != 1
    or receipt.get("logical_tasks") != 22
    or receipt.get("requested_gpu_allocations") != 1
    or receipt.get("central_queue_registration_present") is not True
    or receipt.get("central_queue_unclaimed_during_validation") is not True
    or receipt.get("admission_gate_present_during_validation") is not False
    or receipt.get("direct_gpu_submission_authorized") is not False
    or receipt.get("partial_ranking_authorized") is not False
    or receipt.get("development_evaluator_open") is not False
    or receipt.get("prediction_values_read") is not False
    or receipt.get("metric_values_read") is not False
    or receipt.get("raw_outcomes_read") is not False
    or receipt.get("sealed_labels_or_outcomes_read") is not False
    or receipt.get("gpu_job_submitted_by_admission") is not False
    or preflight.get("status") != "pass"
    or preflight.get("bundle_id") != "model-training-604"
    or preflight.get("fits") != 11
    or preflight.get("logical_tasks") != 22
    or preflight.get("direct_gpu_submission_performed") is not False
    or len(fits) != 11
    or {(int(row["outer_fold"]), int(row["seed"])) for row in fits}
    != {
        (0, 20260824), (0, 20260825), (0, 20260826),
        (1, 20260825), (1, 20260826),
        (2, 20260825), (2, 20260826),
        (3, 20260825), (3, 20260826),
        (4, 20260825), (4, 20260826),
    }
    or "Ran 5 tests" not in stderr
    or not stderr.rstrip().endswith("OK")
):
    raise SystemExit("scBasset three-seed admission receipt differs")

freeze_tree(
    root,
    {
        "artifact_class": "scbasset_three_seed_rectangle_admission",
        "model_id": "scbasset",
        "dataset_id": "gse296875",
        "slurm_job_id": job_id,
        "tests": 5,
        "expected_fits": 15,
        "complete_fits_at_freeze": 4,
        "missing_fits_admitted": 11,
        "central_bundle_id": "model-training-604",
        "central_bundle_count": 1,
        "logical_tasks": 22,
        "requested_gpu_allocations": 1,
        "requested_gpu_hours": 8,
        "queue_gate_materialized_by_this_artifact": True,
        "direct_gpu_submission_performed": False,
        "partial_ranking_performed": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "raw_outcomes_read": False,
        "sealed_labels_or_outcomes_read": False,
        "contract_sha256": receipt["contract_sha256"],
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'scBasset three-seed rectangle admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
