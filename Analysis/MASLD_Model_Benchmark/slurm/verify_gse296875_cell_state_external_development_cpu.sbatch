#!/usr/bin/env bash
#SBATCH --job-name=model-check-331
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT}/Analysis/MASLD_Model_Benchmark
SOURCE_ROOT=${PROJECT}/Analysis/SingleCell/continual_integration
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
RUNTIME_PREFIX=/gpfs/commons/home/jameslee/micromamba/envs/immune-atlas
RUNTIME_LOCK=${ROOT}/executions/environment-locks/cpu_baseline_smoke--8a7dfd2a15a635d5/environment.lock.json
RUNTIME_LOCK_SHA256=8a7dfd2a15a635d5ef5b25fa8ee68355e6c6cfd047ec464cdc1dea83b7201298
SOURCE_LOCK=${SOURCE_ROOT}/reference/gse296875_source_lock_v44.json
MAPPING=${SOURCE_ROOT}/candidates/external-gse296875-metadata-v45-2026-08-17T230000Z/mapping_metadata.tsv
LABELS=${SOURCE_ROOT}/candidates/external-gse296875-author-labels-v46-2026-08-17T234000Z/author_labels.tsv
SCRIPT=${ROOT}/scripts/verify_gse296875_cell_state_external_development.py
TEST=${ROOT}/tests/unit/test_verify_gse296875_cell_state_external_development.py
WRAPPER=${ROOT}/slurm/verify_gse296875_cell_state_external_development_cpu.sbatch
CONFIG=${ROOT}/config/campaigns/gse296875_cell_state_external_development_verification_20260825.json
TARGET=${ROOT}/executions/model-check-331-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite independent verification artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src

module purge
module load python/3.11.5-GCCcore-13.2.0
eval "$(python - "${CONFIG}" <<'PY'
from pathlib import Path
import json
import shlex
import sys
value = json.loads(Path(sys.argv[1]).read_text())
print('PRODUCER=' + shlex.quote(value['producer']['path']))
print('PRODUCER_SHA256=' + shlex.quote(value['producer']['artifacts_sha256']))
PY
)"
python - "${CONFIG}" "${WRAPPER}" "${SCRIPT}" "${TEST}" <<'PY'
from pathlib import Path
import hashlib
import json
import sys
config = json.loads(Path(sys.argv[1]).read_text())
if config.get("schema_version") != "masld-bench-gse296875-cell-state-external-development-verification-campaign-v1":
    raise SystemExit("verification campaign schema differs")
for path in map(Path, sys.argv[2:]):
    expected = config["source_bindings"].get(path.name)
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != expected:
        raise SystemExit(f"verification source differs: {path}")
PY
printf '%s  %s\n' "${PRODUCER_SHA256}" "${PRODUCER}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_LOCK_SHA256}" "${RUNTIME_LOCK}" | sha256sum --check --strict
python - "${PRODUCER}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root)
verify_frozen_tree(root / 'view')
verify_frozen_tree(root / 'view/features')
verify_frozen_tree(root / 'view/evaluator')
PY
module purge

"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH=${ROOT}:${ROOT}/src \
  python -m unittest -v tests.unit.test_verify_gse296875_cell_state_external_development \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH=${ROOT}:${ROOT}/src \
  python "${SCRIPT}" \
    --producer "${PRODUCER}" \
    --expected-producer-sha256 "${PRODUCER_SHA256}" \
    --source-lock "${SOURCE_LOCK}" \
    --mapping-metadata "${MAPPING}" \
    --author-labels "${LABELS}" \
    --output "${STAGE}/verification" \
    > "${STAGE}/validation/run_stdout.json" \
    2> "${STAGE}/validation/run_stderr.txt"
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env PYTHONNOUSERSITE=1 \
  python -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"

module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}" "${PRODUCER_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
import hashlib
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / 'verification/verification_receipt.json').read_text())
if (
    receipt.get('status') != 'pass_external_development_view_only'
    or receipt.get('raw_counts_rederived_exactly') is not True
    or receipt.get('selection_rederived_independently') is not True
    or receipt.get('rows') != 7500
    or receipt.get('donors') != 39
    or receipt.get('external_or_champion_claim_allowed') is not False
    or receipt.get('sealed_outcomes_read') is not False
):
    raise SystemExit('independent verification receipt differs')
freeze_tree(root / 'verification', {
    'artifact_class': 'gse296875_rna_cell_state_external_development_verification',
    'producer_artifacts_sha256': sys.argv[2],
    'rows': 7500,
    'donors': 39,
    'external_or_champion_claim_allowed': False,
    'sealed_outcomes_accessed': False,
    'status': 'passed',
})
verify_frozen_tree(root / 'verification')
verification_sha256 = hashlib.sha256(
    (root / 'verification/ARTIFACTS.json').read_bytes()
).hexdigest()
freeze_tree(root, {
    'artifact_class': 'gse296875_rna_cell_state_external_development_verification_campaign',
    'producer_artifacts_sha256': sys.argv[2],
    'verification_artifacts_sha256': verification_sha256,
    'sealed_outcomes_accessed': False,
    'status': 'passed',
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/verification/ARTIFACTS.json"
