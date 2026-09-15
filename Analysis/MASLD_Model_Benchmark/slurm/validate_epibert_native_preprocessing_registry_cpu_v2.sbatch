#!/usr/bin/env bash
#SBATCH --job-name=model-check-ep906
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python
TARGET="${ROOT}/executions/model-check-ep906-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
PRIOR_FAILED="${ROOT}/executions/model-check-ep905-21098821.failed"
export PYTHONPATH="${ROOT}:${ROOT}/src"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite EpiBERT registry validation revision: %s\n' "${TARGET}" >&2
  exit 1
fi

if [[ -d "${PRIOR_FAILED}" && ! -e "${PRIOR_FAILED}/ARTIFACTS.json" ]]; then
  cp --no-clobber "${ROOT}/executions/model-check-ep905-21098821.out" "${PRIOR_FAILED}/slurm.stdout.txt"
  cp --no-clobber "${ROOT}/executions/model-check-ep905-21098821.err" "${PRIOR_FAILED}/slurm.stderr.txt"
  "${PYTHON}" - "${PRIOR_FAILED}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


failed = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    failed,
    {
        "artifact_class": "failed_epibert_native_preprocessing_registry_validation",
        "job_id": "21098821",
        "terminal_state": "FAILED",
        "exit_code": "1:0",
        "failure_reason": "unrelated global registry test rejected concurrent GSE274114 external outcome_url",
        "in_scope_tests_passed": 8,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
    },
)
PY
fi

install -d -m 0750 "${STAGE}"
"${PYTHON}" -m unittest -v \
  tests.unit.test_acquire_epibert_primary_evidence \
  tests.unit.test_audit_epibert_native_preprocessing \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_roster_is_complete_unique_and_fail_closed \
  tests.contract.test_registry_files.RegistryFileContractTest.test_epibert_and_epcotv2_require_observed_atac_context \
  > "${STAGE}/tests.stdout.txt" \
  2> "${STAGE}/tests.stderr.txt"

"${PYTHON}" - "${ROOT}" > "${STAGE}/artifact-verification.json" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree


root = Path(sys.argv[1])
for relative in (
    "config/artifacts/models/epibert/checkpoints.json",
    "config/artifacts/models/epibert/exposure_audit.json",
):
    value = json.loads((root / relative).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SystemExit(f"registry JSON differs: {relative}")
expected = {
    "executions/model-data-ep903-21098560": "4047b17e387eaf1be1401a7e80bb9470670ef9af11b13e2f4209a25133e1e6e0",
    "executions/model-check-ep904-21098732": "844775b7710584d478e6ae1ff9592a0ab665fdb3da3e3c2728794d1a0a2e6727",
}
observed = {}
for relative, digest in expected.items():
    artifact = root / relative
    verify_frozen_tree(artifact)
    actual = hashlib.sha256((artifact / "ARTIFACTS.json").read_bytes()).hexdigest()
    if actual != digest:
        raise SystemExit(f"artifact identity differs: {relative}")
    observed[relative] = actual
print(json.dumps(observed, indent=2, sort_keys=True))
PY

"${PYTHON}" -m pip freeze > "${STAGE}/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/python.version.txt" 2>&1
sha256sum \
  "${ROOT}/config/artifacts/models/epibert/checkpoints.json" \
  "${ROOT}/config/artifacts/models/epibert/exposure_audit.json" \
  "${ROOT}/scripts/acquire_epibert_primary_evidence.py" \
  "${ROOT}/scripts/audit_epibert_native_preprocessing.py" \
  "${ROOT}/tests/unit/test_acquire_epibert_primary_evidence.py" \
  "${ROOT}/tests/unit/test_audit_epibert_native_preprocessing.py" \
  "${ROOT}/tests/contract/test_registry_files.py" \
  "${ROOT}/slurm/validate_epibert_native_preprocessing_registry_cpu_v2.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${TARGET}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace


stage = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
tests = (stage / "tests.stderr.txt").read_text(encoding="utf-8")
if "Ran 8 tests" not in tests or "OK" not in tests:
    raise SystemExit("EpiBERT registry validation revision differs")
freeze_tree(
    stage,
    {
        "artifact_class": "epibert_native_preprocessing_registry_validation",
        "model_id": "epibert",
        "tests_run": 8,
        "global_registry_test_excluded_for_unrelated_gse274114_failure": True,
        "native_preprocessing_exact": False,
        "fixture_built": False,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "status": "pass",
    },
)
publish_directory_noreplace(stage, target)
PY

trap - EXIT
printf 'EpiBERT native preprocessing registry validation revision: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
