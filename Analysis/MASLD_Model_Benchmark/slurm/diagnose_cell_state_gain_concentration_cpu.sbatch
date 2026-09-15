#!/usr/bin/env bash
#SBATCH --job-name=model-check-902
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCORE=${ROOT}/executions/model-cpu-score-502-21097344/evaluation
SCORE_SHA256=$(sha256sum "${SCORE}/ARTIFACTS.json" | cut -d' ' -f1)
TARGET=${ROOT}/executions/model-check-902-${SLURM_JOB_ID}-gain-concentration
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite diagnostic artifact: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}/diagnostic" "${STAGE}/validation"
cd "${ROOT}"
hostname; "${PYTHON}" -V

# Prove the read-only input is unchanged, before and after.
sha256sum "${SCORE}/ARTIFACTS.json" "${SCORE}/study_summaries.tsv" \
  > "${STAGE}/validation/inputs_before.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_stats tests.unit.test_auprc_reference \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" scripts/diagnose_cell_state_gain_concentration.py \
  --study-summaries "${SCORE}/study_summaries.tsv" \
  --output "${STAGE}/diagnostic/gain_concentration.json" \
  > "${STAGE}/diagnostic/stdout.json"

sha256sum "${SCORE}/ARTIFACTS.json" "${SCORE}/study_summaries.tsv" \
  > "${STAGE}/validation/inputs_after.sha256"
diff "${STAGE}/validation/inputs_before.sha256" "${STAGE}/validation/inputs_after.sha256"
printf 'read-only inputs byte-identical before and after\n'

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum \
  scripts/diagnose_cell_state_gain_concentration.py \
  src/masld_bench/evaluators/stats.py \
  src/masld_bench/evaluators/auprc_reference.py \
  slurm/diagnose_cell_state_gain_concentration_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
payload = json.loads((root / "diagnostic" / "gain_concentration.json").read_text())
assert payload["is_post_hoc"] is True, "a post-hoc diagnostic must say so"
assert payload["revises_nothing"] is True
assert payload["supports_pass_fail_verdict"] is False
freeze_tree(root, {
    "artifact_class": "cell_state_gain_concentration_post_hoc_diagnostic",
    "is_post_hoc": True,
    "revises_nothing": True,
    "supports_pass_fail_verdict": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'diagnostic: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
