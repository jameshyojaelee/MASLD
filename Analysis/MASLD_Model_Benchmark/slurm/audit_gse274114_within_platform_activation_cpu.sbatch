#!/usr/bin/env bash
#SBATCH --job-name=data-check-238
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-check-238-%j.out
#SBATCH --error=executions/logs/data-check-238-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/data-check-238-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE274114 within-platform audit.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

"${PYTHON}" -m unittest -v \
  tests.unit.test_audit_gse274114_within_platform_activation_unexposed \
  > "${STAGE}/unittest.stdout.txt" 2> "${STAGE}/unittest.stderr.txt"
"${PYTHON}" - <<'PY'
from pathlib import Path

for relative in (
    "scripts/audit_gse274114_within_platform_activation_unexposed.py",
    "tests/unit/test_audit_gse274114_within_platform_activation_unexposed.py",
):
    path = Path(relative)
    compile(path.read_text(encoding="utf-8"), path.as_posix(), "exec")
PY
bash -n slurm/audit_gse274114_within_platform_activation_cpu.sbatch
"${PYTHON}" scripts/audit_gse274114_within_platform_activation_unexposed.py \
  --fixture-root "${ROOT}/executions/model-data-080-21082249-gse274114" \
  --dataset "${ROOT}/config/datasets/gse274114_mash_hbv.toml" \
  --hiseq-task "${ROOT}/config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml" \
  --novaseq-task "${ROOT}/config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml" \
  --promotion-gate "${ROOT}/config/evaluation/gse274114_within_instrument_etiology_promotion_gate.json" \
  --output "${STAGE}/audit" \
  > "${STAGE}/audit.stdout.json" 2> "${STAGE}/audit.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/audit_gse274114_within_platform_activation_unexposed.py \
  tests/unit/test_audit_gse274114_within_platform_activation_unexposed.py \
  slurm/audit_gse274114_within_platform_activation_cpu.sbatch \
  config/datasets/gse274114_mash_hbv.toml \
  config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml \
  config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml \
  config/evaluation/gse274114_within_instrument_etiology_promotion_gate.json \
  executions/model-data-080-21082249-gse274114/ARTIFACTS.json \
  "${STAGE}/audit/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
verify_frozen_tree(root / "audit")
receipt = json.loads((root / "audit/receipt.json").read_text(encoding="utf-8"))
tests = (root / "unittest.stderr.txt").read_text(encoding="utf-8")
if (
    "Ran 4 tests" not in tests
    or "OK" not in tests
    or receipt.get("status") != "passed_outcome_blind_within_platform_task_activation"
    or receipt.get("evaluator_label_files_opened") is not False
    or receipt.get("evaluator_fold_file_opened") is not False
    or receipt.get("molecular_quant_values_opened") is not False
    or receipt.get("participant_label_map_read") is not False
    or receipt.get("model_fit_or_scoring_run") is not False
    or receipt.get("global_four_class_task_activated") is not False
    or receipt.get("global_ood_task_activated") is not False
    or receipt.get("sealed_or_champion_claim_eligible") is not False
):
    raise SystemExit("GSE274114 outcome-blind activation firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse274114_within_platform_activation_audit_execution",
        "status": receipt["status"],
        "evaluator_labels_opened": False,
        "molecular_values_opened": False,
        "model_fit_or_scoring_run": False,
        "global_ood_task_activated": False,
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/audit/ARTIFACTS.json"
