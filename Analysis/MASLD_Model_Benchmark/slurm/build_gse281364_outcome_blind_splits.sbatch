#!/usr/bin/env bash
#SBATCH --job-name=bench-mpra-splits
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/gse281364-validated-reconstruction-21066470"
SOURCE_SHA=2a37531b9909201a9b42ff817ca92ff2281b58672ffa6bd36d7dc2e86d3765ba
TARGET="${ROOT}/executions/gse281364-outcome-blind-splits-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite MPRA split artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' "${SOURCE_SHA}" "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
PY
install -d -m 0750 "${STAGE}/validation"
PYTHONPATH="${ROOT}" python -m unittest -v \
  tests.unit.test_build_gse281364_outcome_blind_splits \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="${ROOT}" python "${ROOT}/scripts/build_gse281364_outcome_blind_splits.py" \
  --elements "${SOURCE}/validated/elements.tsv" \
  --source-artifacts "${SOURCE}/ARTIFACTS.json" \
  --source-artifacts-sha256 "${SOURCE_SHA}" \
  --output "${STAGE}/split" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/build_gse281364_outcome_blind_splits.py" \
  "${ROOT}/tests/unit/test_build_gse281364_outcome_blind_splits.py" \
  "${ROOT}/slurm/build_gse281364_outcome_blind_splits.sbatch" \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "split/receipt.json").read_text())
if (
    receipt["status"] != "pass_outcome_blind_split_contract"
    or receipt["qualified_paired_snv_elements"] != 4359
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["outer_folds"] != 5
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
):
    raise SystemExit("GSE281364 outcome-blind split receipt differs")
freeze_tree(root, {
    "artifact_class": "gse281364_outcome_blind_variant_split",
    "dataset_id": "gse281364",
    "qualified_paired_snv_elements": 4359,
    "outer_locus_sequence_groups": 1033,
    "outer_folds": 5,
    "outcomes_read": False,
    "champion_eligible": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 outcome-blind split: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
