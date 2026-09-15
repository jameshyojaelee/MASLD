#!/usr/bin/env bash
#SBATCH --job-name=data-check-251
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/logs/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/logs/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/data-check-251-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite external-cohort readiness validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
bash -n slurm/validate_external_cohort_development_readiness_cpu.sbatch
bash -n slurm/freeze_gse105127_development_readiness_cpu.sbatch
bash -n slurm/freeze_gse268273_development_readiness_cpu.sbatch
"${PYTHON}" -m unittest \
  tests.unit.test_freeze_external_cohort_development_readiness \
  > "${STAGE}/test_stdout.txt" 2> "${STAGE}/test_stderr.txt"
"${PYTHON}" -m py_compile \
  scripts/freeze_external_cohort_development_readiness.py \
  tests/unit/test_freeze_external_cohort_development_readiness.py
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/freeze_external_cohort_development_readiness.py \
  tests/unit/test_freeze_external_cohort_development_readiness.py \
  slurm/validate_external_cohort_development_readiness_cpu.sbatch \
  slurm/freeze_gse105127_development_readiness_cpu.sbatch \
  slurm/freeze_gse268273_development_readiness_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "external_cohort_development_readiness_validation",
        "cohorts": ["GSE105127", "GSE268273"],
        "outcomes_read": False,
        "fit_or_score_performed": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
