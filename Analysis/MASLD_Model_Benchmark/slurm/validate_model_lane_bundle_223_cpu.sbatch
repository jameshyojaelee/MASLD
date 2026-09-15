#!/usr/bin/env bash
#SBATCH --job-name=model-check-223
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=executions/logs/model-check-223-%j.out
#SBATCH --error=executions/logs/model-check-223-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/model-check-223-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
cd "${ROOT}"
if [[ -e ${TARGET} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite model-lane validation bundle.\n' >&2
  exit 1
fi
mkdir -p "${STAGE}/validation"
(
  bash slurm/validate_model_lane_202_cpu.sbatch
) >"${STAGE}/validation/contracts.stdout.txt" \
  2>"${STAGE}/validation/contracts.stderr.txt"
(
  bash slurm/audit_model_lane_202_cpu.sbatch
) >"${STAGE}/validation/audit.stdout.txt" \
  2>"${STAGE}/validation/audit.stderr.txt"
printf 'passed\tcontracts\npassed\taudit\n' > "${STAGE}/validation/checks.tsv"
sha256sum \
  slurm/validate_model_lane_bundle_223_cpu.sbatch \
  slurm/validate_model_lane_202_cpu.sbatch \
  slurm/audit_model_lane_202_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" \
  "${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python" \
  - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1])
checks = (stage / "validation/checks.tsv").read_text(encoding="utf-8").splitlines()
if checks != ["passed\tcontracts", "passed\taudit"]:
    raise SystemExit("model-lane validation bundle is incomplete")
freeze_tree(stage, {
    "artifact_class": "sequence_model_lane_validation_bundle",
    "models": ["enformer_crested_restricted_port", "sei"],
    "new_gpu_jobs_created": 0,
    "outcomes_read": False,
    "sealed_labels_read": False,
    "status": "passed",
})
verify_frozen_tree(stage)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
