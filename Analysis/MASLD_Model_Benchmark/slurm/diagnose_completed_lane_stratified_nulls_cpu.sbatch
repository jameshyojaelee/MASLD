#!/usr/bin/env bash
#SBATCH --job-name=model-check-905
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
TARGET=${ROOT}/executions/model-check-905-${SLURM_JOB_ID}-stratified-nulls
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}/diagnostic" "${STAGE}/validation"
cd "${ROOT}"
hostname; "${PYTHON}" -V

INPUTS=(
  executions/model-predict-865-21130175/predictions/ARTIFACTS.json
  executions/model-eval-866-21130184/evaluator_only/ARTIFACTS.json
  executions/model-predict-886-21130315-nash_vs_nafl/predictions/ARTIFACTS.json
  executions/model-data-888-21130317-both_arms-labels/evaluator_only/ARTIFACTS.json
  executions/model-data-072-21080651-participants/ARTIFACTS.json
  executions/model-data-086-21110619-summarization/ARTIFACTS.json
)
sha256sum "${INPUTS[@]}" > "${STAGE}/validation/inputs_before.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_auprc_reference tests.unit.test_stats \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" \
  scripts/diagnose_completed_lane_stratified_nulls.py \
  --root "${ROOT}" \
  --output "${STAGE}/diagnostic/stratified_nulls.json" \
  | tee "${STAGE}/diagnostic/stdout.txt"

sha256sum "${INPUTS[@]}" > "${STAGE}/validation/inputs_after.sha256"
diff "${STAGE}/validation/inputs_before.sha256" "${STAGE}/validation/inputs_after.sha256"
printf 'every read input byte-identical before and after\n'

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum \
  scripts/diagnose_completed_lane_stratified_nulls.py \
  src/masld_bench/evaluators/auprc_reference.py \
  slurm/diagnose_completed_lane_stratified_nulls_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
payload = json.loads((root / "diagnostic" / "stratified_nulls.json").read_text())
assert payload["is_post_hoc"] is True and payload["revises_nothing"] is True
for lane in payload["lanes"]:
    assert lane["stratum"] == "scan_date"
    assert lane["models"], "a lane reported no models"
freeze_tree(root, {
    "artifact_class": "completed_lane_stratified_null_post_hoc_diagnostic",
    "is_post_hoc": True,
    "revises_nothing": True,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'diagnostic: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
