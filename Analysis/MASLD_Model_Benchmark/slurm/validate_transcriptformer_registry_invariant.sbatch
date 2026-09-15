#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-contract-258
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/transcriptformer-registry-invariant-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer registry validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
set +e
PYTHONPATH=src python -m unittest -v \
  tests.contract.test_registry_files.RegistryFileContractTest.test_transcriptformer_preflight_fails_closed_on_exposure \
  > "${STAGE}/unittest.stdout.txt" \
  2> "${STAGE}/unittest.stderr.txt"
TEST_EXIT=$?
set -e
if [[ "${TEST_EXIT}" -ne 0 ]]; then
  printf 'Targeted registry/campaign tests failed with exit %s.\n' "${TEST_EXIT}" >&2
  exit "${TEST_EXIT}"
fi
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import tomllib
import sys

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.contracts import ModelManifest
from masld_bench.hashing import sha256_file

root = Path.cwd()
stage = Path(sys.argv[1]).resolve(strict=True)
family = tomllib.loads((root / "config/models/cell_foundation.toml").read_text())
entry = next(item for item in family["models"] if item["model_id"] == "transcriptformer_tf_sapiens")
manifest = ModelManifest.from_family_entry(family, entry)
if manifest.exposure_status.value != "unknown" or not manifest.admission_blocking:
    raise SystemExit("TranscriptFormer unknown-exposure invariant differs")
write_json_exclusive(stage / "validation_receipt.json", {
    "schema_version": "masld-bench-transcriptformer-registry-invariant-validation-v1",
    "model_id": manifest.model_id,
    "exposure_status": manifest.exposure_status.value,
    "admission_blocking": manifest.admission_blocking,
    "registry_sha256": sha256_file(root / "config/models/cell_foundation.toml"),
    "registry_test_sha256": sha256_file(root / "tests/contract/test_registry_files.py"),
    "targeted_test_modules": [
        "tests.contract.test_registry_files.RegistryFileContractTest.test_transcriptformer_preflight_fails_closed_on_exposure",
    ],
    "outcomes_opened": False,
    "status": "pass",
})
freeze_tree(stage, {
    "artifact_class": "transcriptformer_registry_invariant_validation",
    "model_id": manifest.model_id,
    "exposure_status": manifest.exposure_status.value,
    "admission_blocking": manifest.admission_blocking,
    "outcomes_opened": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
