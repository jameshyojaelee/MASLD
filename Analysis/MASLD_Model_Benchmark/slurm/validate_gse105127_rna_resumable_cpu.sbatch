#!/usr/bin/env bash
#SBATCH --job-name=data-check-221
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=executions/logs/data-check-221-%j.out
#SBATCH --error=executions/logs/data-check-221-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
TARGET=${ROOT}/executions/data-check-221-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
cd "${ROOT}"
if [[ -e ${TARGET} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite resumable RNA validation.\n' >&2
  exit 1
fi
mkdir -p "${STAGE}/validation"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
"${RUNTIME}/env/bin/python" -m unittest -v \
  tests.unit.test_download_gse105127_rna_bundle_resumable \
  >"${STAGE}/validation/tests.stdout.txt" \
  2>"${STAGE}/validation/tests.stderr.txt"
"${RUNTIME}/env/bin/python" - <<'PY' >"${STAGE}/validation/remote_identity.txt"
from pathlib import Path
from scripts.download_gse105127_rna_bundle_resumable import validate_remote_identity
validate_remote_identity(Path("executions/data-check-218-21087748"))
print("remote identity admitted")
PY
"${RUNTIME}/env/bin/python" scripts/download_gse105127_rna_bundle_resumable.py --help \
  >"${STAGE}/validation/help.txt"
"${RUNTIME}/env/bin/python" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1])
freeze_tree(stage, {
    "artifact_class": "gse105127_resumable_rna_validation",
    "network_accessed": False,
    "labels_accessed": False,
    "status": "passed",
})
verify_frozen_tree(stage)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
