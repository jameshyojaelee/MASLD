#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-context-map
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${OUTER_FOLD:?OUTER_FOLD is required}"
: "${CONTEXT_JOB_ID:?CONTEXT_JOB_ID is required}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
CONTEXTS="${ROOT}/executions/corgi-gse296875-context-counts-${CONTEXT_JOB_ID}"
TARGET="${ROOT}/executions/corgi-gse296875-mapper-fold${OUTER_FOLD}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi mapper artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${CONTEXTS}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
PY
"${PYTHON}" -m unittest -v tests.unit.test_fit_corgi_masked_context_mapper
"${PYTHON}" scripts/fit_corgi_masked_context_mapper.py --contexts "${CONTEXTS}" --outer-fold "${OUTER_FOLD}" --output "${STAGE}"
"${PYTHON}" -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum scripts/fit_corgi_masked_context_mapper.py tests/unit/test_fit_corgi_masked_context_mapper.py slurm/fit_corgi_masked_context_mapper.sbatch > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" "${OUTER_FOLD}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root, fold = Path(sys.argv[1]).resolve(strict=True), int(sys.argv[2])
receipt = json.loads((root / "receipt.json").read_text())
if receipt["outer_fold"] != fold or receipt["held_ATAC_or_other_outcomes_used"]:
    raise SystemExit("Corgi mapper receipt differs")
freeze_tree(root, {"artifact_class": "fold_fitted_masked_Corgi_context_mapper", "model_id": "corgi_regular", "dataset_id": "gse296875", "outer_fold": fold, "valid_fold": receipt["valid_fold"], "outcomes_read": False, "status": "passed"})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi masked context mapper: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
