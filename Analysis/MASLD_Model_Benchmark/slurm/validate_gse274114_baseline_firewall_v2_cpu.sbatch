#!/usr/bin/env bash
#SBATCH --job-name=data-check-245
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-check-245-%j.out
#SBATCH --error=executions/logs/data-check-245-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
DIAGNOSTIC=${ROOT}/executions/data-check-244-21092183
DIAGNOSTIC_SHA256=a5d230f6a54e00f575e79919170f8bb313698363fb0cb8de7deb36065a690bd9
TRACE_SHA256=de468a54e3e7d8dbc15fa04ad00f98273783b84ef8200cfb6d407b3f209ba7f1
TARGET=${ROOT}/executions/data-check-245-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE274114 v2 synthetic validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${DIAGNOSTIC_SHA256}" "${DIAGNOSTIC}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${TRACE_SHA256}" "${DIAGNOSTIC}/audit/ARTIFACTS.json" \
  | sha256sum --check --strict

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse274114_baseline_firewall \
  > "${STAGE}/firewall_unittest.stdout.txt" \
  2> "${STAGE}/firewall_unittest.stderr.txt"
"${PYTHON}" -m unittest -v \
  tests.unit.test_run_gse274114_baseline_task_campaign \
  > "${STAGE}/campaign_unittest.stdout.txt" \
  2> "${STAGE}/campaign_unittest.stderr.txt"
"${PYTHON}" - <<'PY'
from pathlib import Path

for relative in (
    "scripts/gse274114_baseline_firewall.py",
    "scripts/run_gse274114_baseline_task_campaign.py",
    "tests/unit/test_gse274114_baseline_firewall.py",
    "tests/unit/test_run_gse274114_baseline_task_campaign.py",
):
    path = Path(relative)
    compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
bash -n slurm/validate_gse274114_baseline_firewall_v2_cpu.sbatch
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/gse274114_baseline_firewall.py \
  scripts/run_gse274114_baseline_task_campaign.py \
  tests/unit/test_gse274114_baseline_firewall.py \
  tests/unit/test_run_gse274114_baseline_task_campaign.py \
  config/campaigns/gse274114_within_platform_baselines_v2.toml \
  config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml \
  config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml \
  slurm/validate_gse274114_baseline_firewall_v2_cpu.sbatch \
  "${DIAGNOSTIC}/ARTIFACTS.json" \
  "${DIAGNOSTIC}/audit/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1])
firewall_tests = (root / "firewall_unittest.stderr.txt").read_text(encoding="utf-8")
campaign_tests = (root / "campaign_unittest.stderr.txt").read_text(encoding="utf-8")
if (
    "Ran 10 tests" not in firewall_tests
    or "OK" not in firewall_tests
    or "Ran 4 tests" not in campaign_tests
    or "OK" not in campaign_tests
    or "ConvergenceWarning" in firewall_tests
    or "FutureWarning" in firewall_tests
    or "ConvergenceWarning" in campaign_tests
    or "FutureWarning" in campaign_tests
):
    raise SystemExit("GSE274114 v2 synthetic tests differ")
freeze_tree(
    root,
    {
        "artifact_class": "gse274114_baseline_firewall_v2_synthetic_validation",
        "status": "passed_synthetic_only",
        "real_evaluator_labels_opened": False,
        "real_evaluator_folds_opened": False,
        "real_molecular_values_opened": False,
        "real_model_fit_or_scoring_run": False,
        "training_only_ranked_convergence_fallback_tested": True,
        "all_outer_candidates_nonconvergent_fail_closed_tested": True,
        "trainer_evaluator_firewall_tested": True,
        "incomplete_prediction_set_fail_closed_tested": True,
        "fold_local_preprocessing_tested": True,
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
