#!/usr/bin/env bash
#SBATCH --job-name=model-data-065
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PARTICIPANTS=${ROOT}/executions/model-data-063-21079860
PARTICIPANTS_ARTIFACTS_SHA256=b6846c72d85133990bfa3852d939a1383471b07b892534ef292aaea3af569ff9
TARGET=${ROOT}/executions/model-data-065-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite microarray raw audit artifacts.\n' >&2
  exit 1
fi
cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
printf '%s  %s\n' "${PARTICIPANTS_ARTIFACTS_SHA256}" "${PARTICIPANTS}/ARTIFACTS.json" \
  | sha256sum --check --strict
python - "${PARTICIPANTS}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
python -m unittest -v tests.unit.test_acquire_audit_microarray_transfer_raw
python scripts/acquire_audit_microarray_transfer_raw.py \
  --participants "${PARTICIPANTS}" \
  --participants-artifacts-sha256 "${PARTICIPANTS_ARTIFACTS_SHA256}" \
  --output "${STAGE}"
sha256sum \
  scripts/acquire_audit_microarray_transfer_raw.py \
  tests/unit/test_acquire_audit_microarray_transfer_raw.py \
  slurm/acquire_audit_microarray_transfer_raw.sbatch \
  "${PARTICIPANTS}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "raw_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("raw_array_acquisition_complete") is not True
    or audit.get("raw_CEL_integrity_complete") is not True
    or audit.get("normalization_run") is not False
    or audit.get("model_training_activated") is not False
):
    raise SystemExit("microarray raw audit firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "microarray_transfer_raw_CEL_audit",
        "cohort_family_ids": [
            "gse31803_gse49541_fibrosis_array",
            "antwerp_inserm_shared",
        ],
        "CEL_records": {"GSE49541": 72, "GSE83452": 231},
        "raw_array_acquisition_complete": True,
        "sample_axis_exact": True,
        "platform_probe_mapping_complete": False,
        "normalization_run": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed_raw_CEL_integrity",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
