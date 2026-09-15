#!/usr/bin/env bash
#SBATCH --job-name=model-data-063
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/model-data-062-21079818
SOURCE_ARTIFACTS_SHA256=4c1460c573285fa53fe6464f23498e9c27c610c4be14cf3df3cc2944282ae96d
TARGET=${ROOT}/executions/model-data-063-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite microarray participant audit artifacts.\n' >&2
  exit 1
fi
cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
python -m unittest -v tests.unit.test_audit_microarray_transfer_participants
python scripts/audit_microarray_transfer_participants.py \
  --source "${SOURCE}" \
  --source-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
  --output "${STAGE}"
sha256sum \
  scripts/audit_microarray_transfer_participants.py \
  tests/unit/test_audit_microarray_transfer_participants.py \
  slurm/audit_microarray_transfer_participants.sbatch \
  "${SOURCE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "participant_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("participant_join_authoritative_from_deposited_geo_fields") is not True
    or audit.get("raw_array_acquisition_complete") is not False
    or audit.get("model_training_activated") is not False
    or audit.get("sealed_outcomes_read") is not False
):
    raise SystemExit("microarray participant audit firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "microarray_transfer_participant_audit",
        "cohort_family_ids": [
            "gse31803_gse49541_fibrosis_array",
            "antwerp_inserm_shared",
        ],
        "participants": {"GSE49541": 72, "GSE83452": 171},
        "records": {"GSE49541": 72, "GSE83452": 231},
        "participant_join_authoritative": True,
        "raw_array_acquisition_complete": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed_participant_topology_and_source_phenotypes",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
