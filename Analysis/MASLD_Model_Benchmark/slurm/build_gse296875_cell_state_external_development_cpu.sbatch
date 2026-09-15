#!/usr/bin/env bash
#SBATCH --job-name=model-data-330
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
SOURCE_LOCK=${SOURCE_ROOT}/reference/gse296875_source_lock_v44.json
MAPPING=${SOURCE_ROOT}/candidates/external-gse296875-metadata-v45-2026-08-17T230000Z/mapping_metadata.tsv
LABELS=${SOURCE_ROOT}/candidates/external-gse296875-author-labels-v46-2026-08-17T234000Z/author_labels.tsv
ATLAS_ONTOLOGY=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture/broad_liver_five_v1.json
ATLAS_ARTIFACTS=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture/ARTIFACTS.json
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
RUNTIME_PREFIX=/gpfs/commons/home/jameslee/micromamba/envs/immune-atlas
RUNTIME_LOCK=${ROOT}/executions/environment-locks/cpu_baseline_smoke--8a7dfd2a15a635d5/environment.lock.json
RUNTIME_LOCK_SHA256=8a7dfd2a15a635d5ef5b25fa8ee68355e6c6cfd047ec464cdc1dea83b7201298
SCRIPT=${ROOT}/scripts/build_gse296875_cell_state_external_development.py
TEST=${ROOT}/tests/unit/test_build_gse296875_cell_state_external_development.py
SOURCE_BUILDER=${ROOT}/scripts/build_gse296875_rna_atac_smoke.py
ARTIFACT_MODULE=${ROOT}/src/masld_bench/artifacts.py
WRAPPER=${ROOT}/slurm/build_gse296875_cell_state_external_development_cpu.sbatch
CONFIG=${ROOT}/config/campaigns/gse296875_cell_state_external_development_20260825.json
TARGET=${ROOT}/executions/model-data-330-${SLURM_JOB_ID}
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
  printf 'Refusing to overwrite external-development artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src

module purge
module load python/3.11.5-GCCcore-13.2.0
python - "${CONFIG}" "${SCRIPT}" "${TEST}" "${SOURCE_BUILDER}" \
  "${ARTIFACT_MODULE}" "${WRAPPER}" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

config = json.loads(Path(sys.argv[1]).read_text())
if config.get("schema_version") != "masld-bench-gse296875-cell-state-external-development-campaign-v1":
    raise SystemExit("campaign schema differs")
for path in map(Path, sys.argv[2:]):
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    expected = config["source_bindings"][path.name]
    if observed != expected:
        raise SystemExit(f"source binding differs: {path}")
for item in config["input_bindings"]:
    path = Path(item["path"])
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != item["sha256"]:
        raise SystemExit(f"input binding differs: {path}")
PY

python scripts/freeze_python_environment.py \
  --micromamba "${MICROMAMBA}" --prefix "${RUNTIME_PREFIX}" \
  --runtime-id cpu_baseline_smoke --output-root "${STAGE}/runtime_probe" \
  > "${STAGE}/validation/runtime_probe.json"
observed_lock=$(find "${STAGE}/runtime_probe" -mindepth 2 -maxdepth 2 \
  -name environment.lock.json -print -quit)
if [[ -z ${observed_lock} ]]; then
  printf 'Runtime probe did not produce an environment lock.\n' >&2
  exit 1
fi
printf '%s  %s\n' "${RUNTIME_LOCK_SHA256}" "${observed_lock}" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_LOCK_SHA256}" "${RUNTIME_LOCK}" \
  | sha256sum --check --strict
module purge

"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH=${ROOT}:${ROOT}/src \
  python -m unittest -v tests.unit.test_build_gse296875_cell_state_external_development \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env \
  PYTHONNOUSERSITE=1 PYTHONPATH=${ROOT}:${ROOT}/src \
  python "${SCRIPT}" \
  --source-lock "${SOURCE_LOCK}" \
  --mapping-metadata "${MAPPING}" \
  --author-labels "${LABELS}" \
  --atlas-ontology "${ATLAS_ONTOLOGY}" \
  --output "${STAGE}/view" \
  > "${STAGE}/validation/run_stdout.json" \
  2> "${STAGE}/validation/run_stderr.txt"
"${MICROMAMBA}" run --prefix "${RUNTIME_PREFIX}" env PYTHONNOUSERSITE=1 \
  python -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}:${ROOT}/src python - "${STAGE}" "${ATLAS_ARTIFACTS}" <<'PY'
from pathlib import Path
import hashlib
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
view = root / "view"
feature = verify_frozen_tree(view / "features")
evaluator = verify_frozen_tree(view / "evaluator")
campaign = verify_frozen_tree(view)
campaign_sha256 = hashlib.sha256((view / "ARTIFACTS.json").read_bytes()).hexdigest()
receipt = json.loads((view / "campaign_receipt.json").read_text())
if (
    feature["metadata"].get("rows") != 7500
    or feature["metadata"].get("labels_present") is not False
    or feature["metadata"].get("donor_ids_present") is not False
    or evaluator["metadata"].get("donors") != 39
    or evaluator["metadata"].get("champion_claim_allowed") is not False
    or campaign["metadata"].get("external_or_champion_claim_allowed") is not False
    or receipt.get("sealed_outcomes_read") is not False
    or receipt.get("model_fit") is not False
    or receipt.get("metric_calculated") is not False
    or receipt.get("donor_fingerprint_overlap_status") != "not_run"
):
    raise SystemExit("external-development receipt differs")
freeze_tree(root, {
    "artifact_class": "gse296875_rna_cell_state_external_development_build",
    "view_artifacts_sha256": campaign_sha256,
    "atlas_source_artifacts_path": sys.argv[2],
    "rows": 7500,
    "donors": 39,
    "sealed_outcomes_accessed": False,
    "status": "passed",
})
verify_frozen_tree(root)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/view/ARTIFACTS.json" \
  "${TARGET}/view/features/ARTIFACTS.json" "${TARGET}/view/evaluator/ARTIFACTS.json"
