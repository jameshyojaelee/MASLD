#!/usr/bin/env bash
#SBATCH --job-name=model-check-222
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=executions/logs/model-check-222-%j.out
#SBATCH --error=executions/logs/model-check-222-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/model-check-222-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
cd "${ROOT}"
if [[ -e ${TARGET} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite sequence admission bundle.\n' >&2
  exit 1
fi
mkdir -p "${STAGE}/validation"

run_check() {
  local label=$1
  local wrapper=$2
  (
    bash "${wrapper}"
  ) >"${STAGE}/validation/${label}.stdout.txt" \
    2>"${STAGE}/validation/${label}.stderr.txt"
  printf 'passed\t%s\t%s\n' "${label}" "${wrapper}" >> "${STAGE}/validation/checks.tsv"
}

run_check borzoi_native_preflight slurm/preflight_gse281364_borzoi_native_execution_cpu.sbatch
run_check alphagenome_comparator slurm/preflight_alphagenome_comparator_cpu.sbatch
run_check corgi_regular_disposition slurm/validate_corgi_regular_smoke_plan_cpu.sbatch

sha256sum \
  slurm/validate_sequence_admission_bundle_222_cpu.sbatch \
  slurm/preflight_gse281364_borzoi_native_execution_cpu.sbatch \
  slurm/preflight_alphagenome_comparator_cpu.sbatch \
  slurm/validate_corgi_regular_smoke_plan_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" \
  "${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python" \
  - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1])
checks = (stage / "validation/checks.tsv").read_text(encoding="utf-8").splitlines()
if len(checks) != 3 or not all(line.startswith("passed\t") for line in checks):
    raise SystemExit("sequence admission bundle is incomplete")
freeze_tree(stage, {
    "artifact_class": "sequence_model_admission_validation_bundle",
    "checks": len(checks),
    "checkpoint_forwards_executed": False,
    "outcomes_read": False,
    "sealed_assets_read": False,
    "status": "passed",
})
verify_frozen_tree(stage)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
