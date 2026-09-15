#!/usr/bin/env bash
#SBATCH --job-name=model-training-400
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
: "${BUNDLE_PLAN_ROOT:?frozen bundle plan is required}"
: "${BUNDLE_PLAN_ARTIFACTS_SHA256:?bundle plan artifact hash is required}"
: "${BUNDLE_ID:?bundle identity is required}"
: "${BUNDLE_REQUESTED_TIME:?bundle time contract is required}"
: "${SLURM_JOB_ID:?SLURM_JOB_ID is required}"
: "${BUNDLE_PREFLIGHT_ONLY:=0}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != sequence_gpu_dispatcher_v2 ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 2
fi
if [[ ! "${BUNDLE_ID}" =~ ^sequence_gpu_bundle_[0-9]{3}(_[A-Za-z0-9][A-Za-z0-9_-]*)?$ ]]; then
  printf 'Invalid generic bundle identity.\n' >&2
  exit 2
fi
if [[ ! "${BUNDLE_REQUESTED_TIME}" =~ ^([0-4][0-9]):00:00$ ]]; then
  printf 'Invalid bundle time contract.\n' >&2
  exit 2
fi
if [[ "${BUNDLE_PREFLIGHT_ONLY}" != 0 && "${BUNDLE_PREFLIGHT_ONLY}" != 1 ]]; then
  printf 'Invalid bundle preflight mode.\n' >&2
  exit 2
fi

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTROLLER="${ROOT}/scripts/run_sequence_gpu_bundle.py"
if [[ "${BUNDLE_PREFLIGHT_ONLY}" == 1 ]]; then
  : "${BUNDLE_PREFLIGHT_STAGE:?CPU preflight stage is required}"
  if [[ "${BUNDLE_PREFLIGHT_STAGE}" != "${ROOT}"/executions/*.staging/validation/* ]]; then
    printf 'Bundle preflight stage is outside the admitted validation namespace.\n' >&2
    exit 2
  fi
  STAGE="${BUNDLE_PREFLIGHT_STAGE}"
  TARGET=not_applicable_preflight_only
else
  TARGET="${ROOT}/executions/${BUNDLE_ID}-${SLURM_JOB_ID}"
  STAGE="${TARGET}.staging"
  if [[ -d "${TARGET}" ]]; then
    printf 'Bundle allocation already terminal: %s\n' "${TARGET}"
    exit 0
  fi
fi
if [[ "${BUNDLE_PLAN_ROOT}" != "${ROOT}"/executions/sequence-gpu-bundle-plan-* ]]; then
  printf 'Bundle plan path is outside the admitted execution namespace.\n' >&2
  exit 2
fi
printf '%s  %s\n' "${BUNDLE_PLAN_ARTIFACTS_SHA256}" \
  "${BUNDLE_PLAN_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' fa31c4cbb5c7ff49746c19d9a3bb49621cdf3205f9e485c17952e3c567ec078f \
  "${CONTROLLER}" | sha256sum --check --strict
install -d -m 0750 "${STAGE}/logs" "${STAGE}/task_receipts"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BUNDLE_PLAN_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
CONTROLLER_ARGS=(
  --root "${ROOT}"
  --plan-root "${BUNDLE_PLAN_ROOT}"
  --bundle-id "${BUNDLE_ID}"
  --requested-time "${BUNDLE_REQUESTED_TIME}"
  --job-id "${SLURM_JOB_ID}"
  --stage "${STAGE}"
)
if [[ "${BUNDLE_PREFLIGHT_ONLY}" == 1 ]]; then
  CONTROLLER_ARGS+=(--preflight-only)
fi
PYTHONPATH="${ROOT}/src" python "${CONTROLLER}" \
  "${CONTROLLER_ARGS[@]}"
if [[ "${BUNDLE_PREFLIGHT_ONLY}" == 1 ]]; then
  module purge
  printf 'Bundle wrapper CPU preflight passed: %s\n' "${BUNDLE_ID}"
  exit 0
fi

sha256sum \
  "${BUNDLE_PLAN_ROOT}/ARTIFACTS.json" \
  "${BUNDLE_PLAN_ROOT}/plan/bundle_allocations.tsv" \
  "${BUNDLE_PLAN_ROOT}/plan/bundle_tasks.tsv" \
  "${CONTROLLER}" \
  "${ROOT}/slurm/run_sequence_gpu_bundle.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" \
  "${BUNDLE_PLAN_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "bundle_summary.json").read_text(encoding="utf-8"))
freeze_tree(root, {
    "artifact_class": "sequence_gpu_bundle_allocation_receipt",
    "bundle_id": summary["bundle_id"],
    "slurm_job_id": summary["slurm_job_id"],
    "logical_tasks": summary["logical_tasks"],
    "terminal_status": summary["terminal_status"],
    "continue_after_isolated_failure": True,
    "checkpoint_granularity": "logical_task",
    "bundle_plan_artifacts_sha256": sys.argv[2],
    "future_qos": "nslab",
})
PY
module purge
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
