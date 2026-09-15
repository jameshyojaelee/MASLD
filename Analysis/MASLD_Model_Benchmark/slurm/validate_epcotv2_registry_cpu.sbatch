#!/usr/bin/env bash
#SBATCH --job-name=model-check-113
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738"
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON="${RUNTIME}/env/bin/python"
TARGET="${ROOT}/executions/model-check-113-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failure() {
  status=$?
  if [[ -d "${STAGE}" ]]; then
    if [[ "${status}" -ne 0 && ! -e "${TARGET}.failed" ]]; then
      mv "${STAGE}" "${TARGET}.failed"
    else
      rm -rf -- "${STAGE}"
    fi
  fi
  return "${status}"
}
trap preserve_failure EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite EPCOTv2 registry validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
export PYTHONPATH="${ROOT}/src"
cd "${ROOT}"
"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_epcotv2_native_contract \
  tests.contract.test_registry_files.RegistryFileContractTest.test_epibert_and_epcotv2_require_observed_atac_context \
  > "${STAGE}/validation/tests.stdout.txt" 2>&1

"${PYTHON}" - "${ROOT}" <<'PY' > "${STAGE}/validation/registry.json"
from hashlib import sha256
import json
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1])
bundle = root / "config/artifacts/models/epcotv2"
checkpoint = json.loads((bundle / "checkpoints.json").read_text())
crosswalk = json.loads((bundle / "development_crosswalk.json").read_text())
exposure = json.loads((bundle / "exposure_audit.json").read_text())
artifact = root / exposure["native_contract_audit_record"]["path"]
verify_frozen_tree(artifact)
manifest_sha = sha256((artifact / "ARTIFACTS.json").read_bytes()).hexdigest()
crosswalk_sha = sha256((bundle / "development_crosswalk.json").read_bytes()).hexdigest()
assert manifest_sha == exposure["native_contract_audit_record"]["artifacts_sha256"]
assert crosswalk_sha == exposure["development_crosswalk_record"]["sha256"]
assert crosswalk["crosswalk_scope"]["task_scope"] == ["observed_multiome", "variant_to_regulation"]
assert checkpoint["task_fit"]["rna_conditioned_atac"].startswith("not_supported")
assert "donor-by-lineage pseudobulk" in checkpoint["task_fit"]["observed_multiome"]
print(json.dumps({
    "status": "pass",
    "native_contract_artifacts_sha256": manifest_sha,
    "development_crosswalk_sha256": crosswalk_sha,
    "project_data_read": False,
    "outcomes_read": False,
    "sealed_labels_read": False,
}, sort_keys=True))
PY

sha256sum \
  "${ROOT}/config/artifacts/models/epcotv2/checkpoints.json" \
  "${ROOT}/config/artifacts/models/epcotv2/development_crosswalk.json" \
  "${ROOT}/config/artifacts/models/epcotv2/exposure_audit.json" \
  "${ROOT}/tests/contract/test_registry_files.py" \
  "${ROOT}/slurm/validate_epcotv2_registry_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "epcotv2_registry_validation",
        "model": "epcotv2",
        "tests": 6,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_labels_read": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EPCOTv2 registry validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
