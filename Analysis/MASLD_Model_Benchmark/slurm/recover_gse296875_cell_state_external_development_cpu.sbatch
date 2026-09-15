#!/usr/bin/env bash
#SBATCH --job-name=model-data-332
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
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/model-data-330-21100798.failed/view
SOURCE_SHA256=e00c421a8ba7cbd3443377b6702b876bee90ccf17d39407562c3c592eedf13d5
FEATURE_SHA256=1cacb52032beb27969c69ae2d00cfdf3d18b333f1458bb4fb830f5225d79dc6c
EVALUATOR_SHA256=d76d571e824e8a76473ce8ecb0ee3279e0b27fff9ca27acb22a4dd34ef2875cc
ADMISSION_INCIDENT=${ROOT}/config/artifacts/incidents/gse296875_cell_state_external_development_admission_21100789_21100793.json
PACKAGING_INCIDENT=${ROOT}/config/artifacts/incidents/gse296875_cell_state_external_development_packaging_21100798.json
CONFIG=${ROOT}/config/campaigns/gse296875_cell_state_external_development_recovery_20260825.json
WRAPPER=${ROOT}/slurm/recover_gse296875_cell_state_external_development_cpu.sbatch
TARGET=${ROOT}/executions/model-data-332-${SLURM_JOB_ID}
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
  printf 'Refusing to overwrite recovered external-development artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/incidents"
cd "${ROOT}"

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH=${ROOT}:${ROOT}/src
python - "${CONFIG}" "${WRAPPER}" "${ADMISSION_INCIDENT}" \
  "${PACKAGING_INCIDENT}" <<'PY'
from pathlib import Path
import hashlib
import json
import sys
config = json.loads(Path(sys.argv[1]).read_text())
if config.get('schema_version') != 'masld-bench-gse296875-cell-state-external-development-recovery-v1':
    raise SystemExit('recovery config schema differs')
for path in map(Path, sys.argv[2:]):
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != config['source_bindings'][path.name]:
        raise SystemExit(f'recovery source differs: {path}')
PY
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${FEATURE_SHA256}" "${SOURCE}/features/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EVALUATOR_SHA256}" "${SOURCE}/evaluator/ARTIFACTS.json" | sha256sum --check --strict
python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
source = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(source)
verify_frozen_tree(source / 'features')
verify_frozen_tree(source / 'evaluator')
PY

cp -a "${SOURCE}" "${STAGE}/view"
cp "${ADMISSION_INCIDENT}" "${STAGE}/incidents/"
cp "${PACKAGING_INCIDENT}" "${STAGE}/incidents/"
python - "${STAGE}" "${SOURCE_SHA256}" "${FEATURE_SHA256}" \
  "${EVALUATOR_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
source_sha, feature_sha, evaluator_sha = sys.argv[2:]
view = root / 'view'
verify_frozen_tree(view)
verify_frozen_tree(view / 'features')
verify_frozen_tree(view / 'evaluator')
admission = json.loads((root / 'incidents/gse296875_cell_state_external_development_admission_21100789_21100793.json').read_text())
packaging = json.loads((root / 'incidents/gse296875_cell_state_external_development_packaging_21100798.json').read_text())
if (
    admission.get('scientific_script_launched') is not False
    or admission.get('biological_content_parsed') is not False
    or admission.get('gse289173_feature_or_outcome_read') is not False
    or packaging.get('nested_view_artifacts_sha256') != source_sha
    or packaging.get('gse296875_atac_read') is not False
    or packaging.get('metric_calculated') is not False
):
    raise SystemExit('incident disposition differs')
receipt = {
    'schema_version': 'masld-bench-gse296875-cell-state-external-development-recovery-receipt-v1',
    'source_job_id': '21100798',
    'source_nested_view_artifacts_sha256': source_sha,
    'recovered_view_artifacts_sha256': source_sha,
    'feature_artifacts_sha256': feature_sha,
    'evaluator_artifacts_sha256': evaluator_sha,
    'byte_identical_recovery': True,
    'biological_values_recomputed': False,
    'gse296875_atac_read': False,
    'gse289173_feature_or_outcome_read': False,
    'model_fit': False,
    'prediction_emitted': False,
    'metric_calculated': False,
    'external_or_champion_claim_allowed': False,
    'status': 'passed_recovery_only',
}
(root / 'recovery_receipt.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
freeze_tree(root / 'incidents', {
    'artifact_class': 'gse296875_cell_state_external_development_incident_bundle',
    'admission_jobs': ['21100789', '21100793'],
    'scientific_script_launched': False,
    'biological_content_parsed': False,
    'sealed_outcomes_accessed': False,
    'status': 'passed',
})
freeze_tree(root, {
    'artifact_class': 'gse296875_rna_cell_state_external_development_recovered_build',
    'view_artifacts_sha256': source_sha,
    'feature_artifacts_sha256': feature_sha,
    'evaluator_artifacts_sha256': evaluator_sha,
    'byte_identical_recovery': True,
    'external_or_champion_claim_allowed': False,
    'sealed_outcomes_accessed': False,
    'status': 'passed',
})
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/view/ARTIFACTS.json" \
  "${TARGET}/incidents/ARTIFACTS.json"
