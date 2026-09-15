#!/usr/bin/env bash
#SBATCH --job-name=model-check-044
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TEST_PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/decima/bin/python
TARGET=${ROOT}/executions/geneformer-v2-316m-preflight-validation-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failure() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${exit_code}"
}
trap preserve_failure EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Geneformer validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"

bash -n \
  slurm/build_geneformer_v2_316m_compat_runtime_cpu.sbatch \
  slurm/preflight_geneformer_v2_316m_one_batch_l40s.sbatch \
  slurm/validate_geneformer_v2_316m_preflight_cpu.sbatch
"${TEST_PYTHON}" -m unittest -v \
  tests.unit.test_preflight_geneformer_v2_316m_common \
  tests.unit.test_geneformer_v2_316m_slurm_contract \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
sha256sum \
  scripts/preflight_geneformer_v2_316m_common.py \
  tests/unit/test_preflight_geneformer_v2_316m_common.py \
  tests/unit/test_geneformer_v2_316m_slurm_contract.py \
  slurm/build_geneformer_v2_316m_compat_runtime_cpu.sbatch \
  slurm/preflight_geneformer_v2_316m_one_batch_l40s.sbatch \
  slurm/validate_geneformer_v2_316m_preflight_cpu.sbatch \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
tests = (stage / "tests.stderr.txt").read_text(encoding="utf-8")
if "Ran 5 tests" not in tests or "OK" not in tests:
    raise SystemExit("Geneformer preflight validation tests differ")
freeze_tree(
    stage,
    {
        "artifact_class": "geneformer_v2_316m_preflight_source_validation",
        "model_id": "geneformer_v2_316m",
        "tests_passed": 5,
        "gpu_job_submitted": False,
        "labels_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
