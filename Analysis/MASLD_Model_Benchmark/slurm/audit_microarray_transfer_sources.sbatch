#!/usr/bin/env bash
#SBATCH --job-name=model-data-062
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/model-data-062-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite microarray source audit artifacts.\n' >&2
  exit 1
fi
cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python -m unittest -v tests.unit.test_audit_microarray_transfer_sources
python scripts/audit_microarray_transfer_sources.py --output "${STAGE}"
sha256sum \
  scripts/audit_microarray_transfer_sources.py \
  tests/unit/test_audit_microarray_transfer_sources.py \
  slurm/audit_microarray_transfer_sources.sbatch \
  > "${STAGE}/source.sha256"
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "source_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("total_geo_sample_records") != 303
    or audit.get("participant_join_authoritative") is not False
    or audit.get("phenotype_mapping_authoritative") is not False
    or audit.get("model_training_activated") is not False
):
    raise SystemExit("microarray source audit firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "microarray_transfer_source_audit",
        "cohort_family_ids": [
            "gse31803_gse49541_fibrosis_array",
            "antwerp_inserm_shared",
        ],
        "official_geo_sample_records": 303,
        "participant_join_authoritative": False,
        "phenotype_mapping_authoritative": False,
        "raw_array_acquisition_complete": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed_source_metadata_only",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
