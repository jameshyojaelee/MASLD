#!/usr/bin/env bash
#SBATCH --job-name=data-check-241
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-check-241-%j.out
#SBATCH --error=executions/logs/data-check-241-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/data-check-241-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE274114 baseline-firewall validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse274114_baseline_firewall \
  > "${STAGE}/unittest.stdout.txt" 2> "${STAGE}/unittest.stderr.txt"
"${PYTHON}" - <<'PY'
from pathlib import Path

for relative in (
    "scripts/gse274114_baseline_firewall.py",
    "tests/unit/test_gse274114_baseline_firewall.py",
):
    path = Path(relative)
    compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
bash -n slurm/validate_gse274114_baseline_firewall_cpu.sbatch
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/gse274114_baseline_firewall.py \
  tests/unit/test_gse274114_baseline_firewall.py \
  config/campaigns/gse274114_within_platform_baselines_v1.toml \
  config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml \
  config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml \
  slurm/validate_gse274114_baseline_firewall_cpu.sbatch \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1])
tests = (root / "unittest.stderr.txt").read_text(encoding="utf-8")
if (
    "Ran 8 tests" not in tests
    or "OK" not in tests
    or "ConvergenceWarning" in tests
    or "FutureWarning" in tests
):
    raise SystemExit("GSE274114 synthetic baseline-firewall tests differ")
freeze_tree(
    root,
    {
        "artifact_class": "gse274114_baseline_firewall_synthetic_validation",
        "status": "passed_synthetic_only",
        "real_evaluator_labels_opened": False,
        "real_evaluator_folds_opened": False,
        "real_molecular_values_opened": False,
        "real_model_fit_or_scoring_run": False,
        "trainer_evaluator_firewall_tested": True,
        "incomplete_prediction_set_fail_closed_tested": True,
        "fold_local_preprocessing_tested": True,
        "convergence_warning_fail_closed_tested": True,
        "molecular_path_allowlist_tested": True,
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
