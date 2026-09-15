#!/usr/bin/env bash
#SBATCH --job-name=model-check-285
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
CONFIG="${ROOT}/config/scooby_epicardioids_full_backbone_synthetic_probe.json"
WRAPPER="${ROOT}/slurm/probe_scooby_full_backbone_synthetic_l40s.sbatch"
TARGET="${ROOT}/executions/model-check-285-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Scooby GPU queue validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v \
  tests.unit.test_probe_scooby_full_backbone_synthetic \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PREFLIGHT_ONLY=1 bash "${WRAPPER}" \
  > "${STAGE}/validation/wrapper-preflight.stdout.txt" \
  2> "${STAGE}/validation/wrapper-preflight.stderr.txt"
sbatch --test-only "${WRAPPER}" \
  > "${STAGE}/validation/sbatch-test-only.txt" 2>&1
sed -n '1,13p' "${WRAPPER}" > "${STAGE}/validation/exact-gpu-header.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/probe_scooby_full_backbone_synthetic.py" \
  "${ROOT}/tests/unit/test_probe_scooby_full_backbone_synthetic.py" \
  "${WRAPPER}" \
  "${ROOT}/slurm/validate_scooby_full_backbone_queue_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" "${WRAPPER}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
wrapper = Path(sys.argv[3]).resolve(strict=True)
tests = (stage / "validation/tests.stderr.txt").read_text(encoding="utf-8")
preflight = (stage / "validation/wrapper-preflight.stdout.txt").read_text(encoding="utf-8")
test_only = (stage / "validation/sbatch-test-only.txt").read_text(encoding="utf-8")
text = wrapper.read_text(encoding="utf-8")
if (
    "Ran 3 tests" not in tests
    or "OK" not in tests
    or "Scooby full-backbone preflight passed" not in preflight
    or not test_only.strip()
    or "#SBATCH --job-name=model-probe-284" not in text
    or "#SBATCH --partition=gpu" not in text
    or "#SBATCH --account=nslab" not in text
    or "#SBATCH --qos=nslab" not in text
    or "#SBATCH --gres=gpu:l40s:1" not in text
    or "#SBATCH --array" in text
    or "--qos=innovation" in text
):
    raise SystemExit("Scooby GPU queue preflight differs")
freeze_tree(
    stage,
    {
        "artifact_class": "restricted_model_full_backbone_gpu_queue_validation",
        "model_id": "scooby_epicardioids",
        "wrapper_job_name": "model-probe-284",
        "partition": "gpu",
        "account": "nslab",
        "qos": "nslab",
        "gpu": "l40s:1",
        "logical_tasks": 1,
        "manual_gpu_submission": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "rna_conditioned_atac_eligible": False,
        "open_champion_eligible": False,
        "status": "pass_ready_for_central_dispatcher",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'Scooby full-backbone GPU queue validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
