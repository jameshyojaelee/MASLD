#!/usr/bin/env bash
#SBATCH --job-name=data-stage-902
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/gse256398-source-h5-21069407"
DONORS="${ROOT}/executions/gse256398-donor-metadata-21069678/donor_metadata.tsv"
TARGET="${ROOT}/executions/gse256398-terminal-disposition-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"

cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a

mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# Both inputs are already-frozen trees; verify before binding either.
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import verify_frozen_tree

for relative in (
    "executions/gse256398-source-h5-21069407",
    "executions/gse256398-donor-metadata-21069678",
):
    verify_frozen_tree(Path("${ROOT}") / relative)
print("GSE256398 frozen source and donor-metadata verification passed")
PY

# Full contract suite, no test excluded. See the PXD051911 job for why the
# check is a no-new-failure baseline rather than a subset of tests.
set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_build_gse256398_terminal_disposition \
  tests.unit.test_build_gse256398_donor_metadata \
  tests.unit.test_check_contract_suite_baseline \
  tests.contract.test_registry_files \
  > "${STAGING}/validation/unittest.txt" 2>&1
set -e
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGING}/validation/unittest.txt" \
  --require-module test_build_gse256398_terminal_disposition \
  --require-module test_registry_files \
  | tee "${STAGING}/validation/contract_baseline.txt"

PYTHONPATH=${ROOT}:${ROOT}/src python scripts/build_gse256398_terminal_disposition.py \
  --source "${SOURCE}" \
  --donor-metadata "${DONORS}" \
  --output "${STAGING}/audit/disposition" \
  > "${STAGING}/validation/build.stdout.txt" \
  2> "${STAGING}/validation/build.stderr.txt"

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true

sha256sum \
  "${ROOT}/config/datasets/gse256398.toml" \
  "${ROOT}/scripts/build_gse256398_terminal_disposition.py" \
  "${ROOT}/tests/unit/test_build_gse256398_terminal_disposition.py" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  "${ROOT}/tests/unit/test_check_contract_suite_baseline.py" \
  "${ROOT}/slurm/build_gse256398_terminal_disposition.sbatch" \
  "${DONORS}" \
  > "${STAGING}/source.sha256"

PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree

digest = freeze_tree(
    Path("${STAGING}"),
    metadata={
        "artifact_class": "gse256398_terminal_disposition",
        "dataset_id": "gse256398",
        "status": "terminal_no_executable_biological_taskspec",
        "unit_of_replication": "donor",
        "sealed_outcomes_read": False,
        "barcode_cell_labels_read": False,
        "models_fit": False,
        "models_scored": False,
        "supervised_accuracy_claim_allowed": False,
        "model_selection_allowed": False,
    },
)
print("ARTIFACTS.json sha256:", digest)
PY

mv "${STAGING}" "${TARGET}"
echo "FROZEN: ${TARGET}"
