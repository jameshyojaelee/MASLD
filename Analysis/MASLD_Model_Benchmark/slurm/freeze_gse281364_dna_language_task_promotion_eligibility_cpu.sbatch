#!/usr/bin/env bash
#SBATCH --job-name=model-check-608
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG=${ROOT}/config/gse281364_dna_language_task_promotion_eligibility.json
TARGET=${ROOT}/executions/model-check-608-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite eligibility audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH=${ROOT}:${ROOT}/src
cd "${ROOT}"

SOURCE_FILES=(
  config/gse281364_dna_language_task_promotion_eligibility.json
  scripts/freeze_gse281364_dna_language_task_promotion_eligibility.py
  tests/unit/test_gse281364_dna_language_task_promotion_eligibility.py
  slurm/freeze_gse281364_dna_language_task_promotion_eligibility_cpu.sbatch
)
sha256sum "${SOURCE_FILES[@]}" > "${STAGE}/source.start.sha256"

"${PYTHON}" -m json.tool "${CONFIG}" > /dev/null
"${PYTHON}" -m unittest -v \
  tests.unit.test_gse281364_dna_language_task_promotion_eligibility \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" scripts/freeze_gse281364_dna_language_task_promotion_eligibility.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/eligibility" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"

"${PYTHON}" - "${STAGE}/eligibility" <<'PY'
import json
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass_performance_blind_task_specific_eligibility"
    or receipt["performance_read"]
    or receipt["prediction_values_read"]
    or receipt["reporter_outcomes_read"]
    or receipt["sealed_assets_read"]
    or receipt["metrics_calculated"]
    or receipt["ranking_created"]
    or receipt["shortlist_created"]
    or receipt["promotion_executed"]
    or manifest["metadata"]["promotion_executed"]
):
    raise SystemExit("eligibility artifact firewall differs")
PY

"${PYTHON}" -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum "${SOURCE_FILES[@]}" > "${STAGE}/source.end.sha256"
cmp "${STAGE}/source.start.sha256" "${STAGE}/source.end.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(root, {
    "artifact_class": "gse281364_dna_language_task_promotion_eligibility_campaign",
    "dataset_id": "gse281364",
    "task_id": "mpra_allelic_direction",
    "models": ["hyenadna", "nucleotide_transformer"],
    "performance_read": False,
    "reporter_outcomes_read": False,
    "sealed_assets_read": False,
    "metrics_calculated": False,
    "promotion_executed": False,
    "status": "pass_performance_blind_task_specific_eligibility",
})
verify_frozen_tree(root)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Frozen task-specific eligibility: %s\n' "${TARGET}"
