#!/usr/bin/env bash
#SBATCH --job-name=model-check-221
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=executions/logs/model-check-221-%j.out
#SBATCH --error=executions/logs/model-check-221-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/model-check-221-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
cd "${ROOT}"
if [[ -e ${TARGET} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite admission validation bundle.\n' >&2
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

run_check gse105127_exact_reference slurm/validate_gse105127_reference_legacy_exact_cpu.sbatch
run_check corgi_plus_disposition slurm/validate_corgi_plus_disposition_cpu.sbatch
run_check scooby_released_validation slurm/validate_model_lane_201_cpu.sbatch
run_check scooby_released_audit slurm/audit_model_lane_201_cpu.sbatch
run_check borzoi_native_preflight slurm/preflight_gse281364_borzoi_native_execution_cpu.sbatch

sha256sum \
  slurm/validate_admission_bundle_220_cpu.sbatch \
  slurm/validate_gse105127_reference_legacy_exact_cpu.sbatch \
  slurm/validate_corgi_plus_disposition_cpu.sbatch \
  slurm/validate_model_lane_201_cpu.sbatch \
  slurm/audit_model_lane_201_cpu.sbatch \
  slurm/preflight_gse281364_borzoi_native_execution_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" \
  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python \
  - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1])
checks = (stage / "validation/checks.tsv").read_text(encoding="utf-8").splitlines()
if len(checks) != 5 or not all(line.startswith("passed\t") for line in checks):
    raise SystemExit("admission validation bundle is incomplete")
freeze_tree(stage, {
    "artifact_class": "model_and_data_admission_validation_bundle",
    "checks": len(checks),
    "sealed_assets_read": False,
    "status": "passed",
})
verify_frozen_tree(stage)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
