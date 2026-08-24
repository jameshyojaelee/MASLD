#!/usr/bin/env bash
#SBATCH --job-name=masld-scooby-cross-review
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/scooby-cross-review-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scooby cross-review artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/scooby_cross_review.py" \
  --project-root "${ROOT}" \
  --output "${STAGE}"
sha256sum \
  "${ROOT}/scripts/scooby_cross_review.py" \
  "${ROOT}/slurm/scooby_run_cross_review.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "scooby_cross_review_receipt.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "scooby_cross_review",
        "models": sorted(receipt["models"]),
        "review_order": receipt["review_order"],
        "checkpoint_forward_allowed": receipt["checkpoint_forward_allowed"],
        "sealed_inference_allowed": receipt["sealed_inference_allowed"],
        "rna_conditioned_atac_models": receipt["rna_conditioned_atac_models"],
        "terminal_disposition": receipt["terminal_disposition"],
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scooby cross-review artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
