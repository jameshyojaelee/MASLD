#!/usr/bin/env bash
#SBATCH --job-name=g256-qc-sens
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/gse256398-qc-membership-21069669
TARGET=${ROOT}/executions/gse256398-qc-sensitivity-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

# A failed post-move check must not leave a directory standing at the
# authoritative name. Frozen trees are the campaign's unit of authority, and a
# tree that failed verification is internally consistent enough to pass
# verify_frozen_tree for whoever binds it next.
PROMOTION_PENDING=0
preserve_failed() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
  if [[ "${PROMOTION_PENDING}" -eq 1 && -e "${TARGET}" && ! -e "${FAILED}" ]]; then
    mv "${TARGET}" "${FAILED}"
  fi
}
trap preserve_failed EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 QC sensitivity: %s\n' "${TARGET}" >&2
  exit 1
fi
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH=${ROOT}/src:${ROOT}
python -m unittest -v tests.unit.test_build_gse256398_qc_sensitivity
python "${ROOT}/scripts/build_gse256398_qc_sensitivity.py" \
  --qc-root "${SOURCE}" \
  --output "${STAGE}"
pip freeze > "${STAGE}/pip-freeze.txt"
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1])
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
freeze_tree(
    root,
    {
        "artifact_class": "gse256398_qc_sensitivity",
        "dataset_id": "gse256398",
        "donors": receipt["donors"],
        "pathological_auto_call_donors": receipt["pathological_auto_call_donors"],
        "sensitivity_retained": receipt["sensitivity_retained"],
        "outcome_independent": True,
        "model_training_activated": False,
    },
)
PY
mv "${STAGE}" "${TARGET}"
PROMOTION_PENDING=1
python - "${TARGET}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("GSE256398 QC sensitivity frozen-tree verification passed")
PY
PROMOTION_PENDING=0
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
