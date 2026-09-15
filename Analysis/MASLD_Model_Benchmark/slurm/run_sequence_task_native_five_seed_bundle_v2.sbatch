#!/usr/bin/env bash
#SBATCH --job-name=model-training-700
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
: "${RECTANGLE_PLAN_ROOT:?frozen rectangle plan is required}"
: "${RECTANGLE_PLAN_ARTIFACTS_SHA256:?frozen rectangle plan hash is required}"
: "${RECTANGLE_BUNDLE_ID:?generic bundle identity is required}"
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != sequence_five_seed_rectangle_v1 ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 2
fi
if [[ ! "${RECTANGLE_BUNDLE_ID}" =~ ^model-training-7[0-9][0-9]$ ]]; then
  printf 'Invalid generic bundle identity.\n' >&2
  exit 2
fi

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTROLLER="${ROOT}/scripts/run_sequence_task_native_five_seed_bundle_v2.py"
if [[ "${RECTANGLE_PLAN_ROOT}" != "${ROOT}"/executions/sequence-five-seed-rectangle-admission-* ]]; then
  printf 'Rectangle plan root is outside the admitted namespace.\n' >&2
  exit 2
fi
printf '%s  %s\n' "${RECTANGLE_PLAN_ARTIFACTS_SHA256}" \
  "${RECTANGLE_PLAN_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' a57f8b222da3e2c6179d9300eb58a1186af91a4b87274281a239efbe7a4af6c7 \
  "${CONTROLLER}" | sha256sum --check --strict

TARGET="${ROOT}/executions/${RECTANGLE_BUNDLE_ID}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -d "${TARGET}" ]]; then
  printf 'Bundle is already terminal: %s\n' "${TARGET}"
  exit 0
fi
if [[ -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite bundle staging path: %s\n' "${STAGE}" >&2
  exit 2
fi
install -d -m 0750 "${STAGE}/logs" "${STAGE}/task_receipts"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python "${CONTROLLER}" \
  --root "${ROOT}" \
  --plan-root "${RECTANGLE_PLAN_ROOT}" \
  --bundle-id "${RECTANGLE_BUNDLE_ID}" \
  --job-id "${SLURM_JOB_ID}" \
  --stage "${STAGE}"
sha256sum \
  "${RECTANGLE_PLAN_ROOT}/ARTIFACTS.json" \
  "${RECTANGLE_PLAN_ROOT}/plan/bundle_allocations.tsv" \
  "${RECTANGLE_PLAN_ROOT}/plan/bundle_tasks.tsv" \
  "${CONTROLLER}" \
  "${ROOT}/slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" \
  "${RECTANGLE_PLAN_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "bundle_summary.json").read_text(encoding="utf-8"))
freeze_tree(root, {
    "artifact_class": "sequence_task_native_five_seed_bundle_receipt",
    "bundle_id": summary["bundle_id"],
    "slurm_job_id": summary["slurm_job_id"],
    "logical_fits": summary["logical_fits"],
    "terminal_status": summary["terminal_status"],
    "checkpoint_granularity": "logical_fit",
    "continue_after_isolated_fit_failure": True,
    "rectangle_plan_artifacts_sha256": sys.argv[2],
    "qos": "nslab",
    "outcome_paths_available_to_fits": False,
    "sealed_paths_available_to_fits": False,
    "prediction_or_metric_values_read": False,
})
PY
module purge
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
