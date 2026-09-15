#!/usr/bin/env bash
#SBATCH --job-name=data-check-217
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=executions/logs/data-check-217-%j.out
#SBATCH --error=executions/logs/data-check-217-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
TARGET=${ROOT}/executions/data-check-217-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
cd "${ROOT}"
if [[ -e ${TARGET} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite exact-reference validation.\n' >&2
  exit 1
fi
mkdir -p "${STAGE}/validation"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_gse105127_reference_legacy_exact \
  >"${STAGE}/validation/tests.stdout.txt" \
  2>"${STAGE}/validation/tests.stderr.txt"
"${RUNTIME}/env/bin/python" scripts/freeze_gse105127_reference_bundle_legacy_exact.py --help \
  >"${STAGE}/validation/help.txt"
"${RUNTIME}/env/bin/python" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1])
freeze_tree(stage, {
    "artifact_class": "gse105127_exact_legacy_reference_validation",
    "status": "passed",
})
verify_frozen_tree(stage)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
