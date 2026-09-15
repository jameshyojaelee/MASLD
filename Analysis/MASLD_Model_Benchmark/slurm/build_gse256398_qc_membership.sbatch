#!/usr/bin/env bash
#SBATCH --job-name=source-qc-256398
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/gse256398-source-h5-21069407
CROSSWALK=${ROOT}/executions/gse256398-gencode49-crosswalk-21069439
PROBE=${ROOT}/executions/gse256398-scrublet-qc-probe-21069497
PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
TARGET=${ROOT}/executions/gse256398-qc-membership-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 QC membership: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${SOURCE}" "${CROSSWALK}" "${PROBE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("GSE256398 source, crosswalk, and QC probe frozen-tree verification passed")
PY
module purge
PYTHONPATH=${ROOT}:${PYTHONPATH:-} "${PYTHON}" -m unittest -v \
  tests.unit.test_gse256398_scrublet_qc
PYTHONPATH=${ROOT}:${PYTHONPATH:-} "${PYTHON}" \
  scripts/build_gse256398_qc_membership.py \
  --source "${SOURCE}" \
  --crosswalk "${CROSSWALK}" \
  --output "${STAGE}" \
  --seed 0
"${PYTHON}" -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum \
  scripts/probe_gse256398_scrublet_qc.py \
  scripts/build_gse256398_qc_membership.py \
  tests/unit/test_gse256398_scrublet_qc.py \
  slurm/build_gse256398_qc_membership.sbatch \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
with (root / "retained_membership.tsv").open(encoding="utf-8", newline="") as handle:
    retained_rows = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
if (
    receipt["status"] != "pass"
    or receipt["donors"] != 26
    or receipt["nuclei_input"] != 197_942
    or receipt["final_retained"] != retained_rows
    or receipt["disease_age_sex_or_histology_read"]
    or receipt["normalization_or_variable_feature_selection_fitted"]
    or receipt["model_training_activated"]
):
    raise SystemExit("GSE256398 full QC membership receipt differs")
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse256398_qc_membership",
            "dataset_id": "gse256398",
            "donors": receipt["donors"],
            "nuclei_input": receipt["nuclei_input"],
            "final_retained": receipt["final_retained"],
            "scrublet_seed": receipt["seed"],
            "barcode_cell_labels": "structurally_missing",
            "outcome_independent": True,
            "model_training_activated": False,
        },
    )
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE256398 full QC membership: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
