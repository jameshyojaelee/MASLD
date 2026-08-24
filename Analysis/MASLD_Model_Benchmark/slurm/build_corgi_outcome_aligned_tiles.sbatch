#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-tiles
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
REFERENCE="${ROOT}/executions/sequence-fasta-21063808"
FASTA="${REFERENCE}/GRCh38.p14.sequence_model.fa"
FAI="${FASTA}.fai"
TARGET="${ROOT}/executions/corgi-outcome-aligned-tiles-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi tile artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f ${SPLIT}/ARTIFACTS.json" \
  "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 ${REFERENCE}/ARTIFACTS.json" \
  "ed20dceeb315cffa2bc3de050d846eceefc552451fef741f0dd95b0149d10fce ${FASTA}" \
  "7c396dc498b1c064852b2e868b7d0fb96327d2dbd5c9fe80f43165bbe9dc0ebb ${FAI}"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/validation"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}" python -m unittest \
  tests.unit.test_build_corgi_outcome_aligned_tiles \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" python "${ROOT}/scripts/build_corgi_outcome_aligned_tiles.py" \
  --split-contract "${SPLIT}" --reference-fasta "${FASTA}" --reference-fai "${FAI}" \
  --output "${STAGE}/tiles" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_corgi_outcome_aligned_tiles.py" \
  "${ROOT}/tests/unit/test_build_corgi_outcome_aligned_tiles.py" \
  "${ROOT}/slurm/build_corgi_outcome_aligned_tiles.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "tiles" / "summary.json").read_text())
if (
    summary["status"] != "pass_outcome_free_tile_contract"
    or summary["outcomes_or_labels_read"]
    or summary["model_predictions_read"]
    or summary["scoreable_windows"] + summary["excluded_windows"] != 80_000
):
    raise SystemExit("Corgi tile receipt differs")
freeze_tree(root, {
    "artifact_class": "corgi_outcome_aligned_tile_contract",
    "dataset_id": "gse296875",
    "tiles": summary["tiles"],
    "scoreable_windows": summary["scoreable_windows"],
    "excluded_windows": summary["excluded_windows"],
    "outcomes_read": False,
    "predictions_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi outcome-aligned tile contract: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
