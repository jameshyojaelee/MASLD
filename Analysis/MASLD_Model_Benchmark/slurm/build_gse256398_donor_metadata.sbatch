#!/usr/bin/env bash
#SBATCH --job-name=source-metadata-256398
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/gse256398-source-h5-21069407
TARGET=${ROOT}/executions/gse256398-donor-metadata-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 donor metadata: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("GSE256398 source frozen-tree verification passed")
PY
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_build_gse256398_donor_metadata
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/build_gse256398_donor_metadata.py \
  --source "${SOURCE}" \
  --output "${STAGE}"
sha256sum \
  scripts/build_gse256398_donor_metadata.py \
  tests/unit/test_build_gse256398_donor_metadata.py \
  slurm/build_gse256398_donor_metadata.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH=${ROOT}/src python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass"
    or receipt["donors"] != 26
    or receipt["age_observed_donors"] != 26
    or receipt["recorded_sex_observed_donors"] != 26
    or receipt["nas_observed_donors"] != 0
    or receipt["sealed_or_final_outcomes_read"]
    or receipt["model_training_activated"]
):
    raise SystemExit("GSE256398 donor metadata receipt differs")
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse256398_donor_metadata",
            "dataset_id": "gse256398",
            "donors": 26,
            "age_observed_donors": 26,
            "recorded_sex_observed_donors": 26,
            "nas_observed_donors": 0,
            "unit_of_replication": "donor",
            "model_training_activated": False,
            "sealed": False,
        },
    )
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE256398 donor metadata: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
