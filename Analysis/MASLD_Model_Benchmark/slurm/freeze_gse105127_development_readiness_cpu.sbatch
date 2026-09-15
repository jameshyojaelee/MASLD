#!/usr/bin/env bash
#SBATCH --job-name=data-check-252
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/logs/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/logs/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${READINESS_VALIDATION_ROOT:?READINESS_VALIDATION_ROOT must name the frozen validation}"
: "${READINESS_VALIDATION_SHA256:?READINESS_VALIDATION_SHA256 must pin that ARTIFACTS.json}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
VALIDATION=${READINESS_VALIDATION_ROOT}
VALIDATION_SHA256=${READINESS_VALIDATION_SHA256}
AUDITOR_SHA256=61d809ee9ac5064ccd443704b063296d233d7b34430f917d3f21d8f9d9183a56
PLAN=${ROOT}/executions/gse105127-production-9675cd48993468d4/plan
PLAN_SHA256=c4f91d2a3fedc7b3b1802a0d6cf3d1b1803ce1bea70cc343f6e18f744783b973
CAMPAIGN=${ROOT}/executions/gse105127-gzip-preflight-v2-4f96d314e84de086
SOURCE_AUDIT=${ROOT}/executions/data-check-225-21089700
SOURCE_AUDIT_SHA256=e24160f09cb985c8be6621b63d2b1389f3ab27c0b8f8601f9eb962e693fa5880
TARGET=${ROOT}/executions/data-check-252-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE105127 readiness disposition.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/disposition"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${VALIDATION_SHA256}" "${VALIDATION}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${AUDITOR_SHA256}" \
  "${ROOT}/scripts/freeze_external_cohort_development_readiness.py" \
  | sha256sum --check --strict
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
bash -n slurm/freeze_gse105127_development_readiness_cpu.sbatch
rmdir "${STAGE}/disposition"
"${PYTHON}" scripts/freeze_external_cohort_development_readiness.py gse105127 \
  --plan-root "${PLAN}" \
  --campaign-root "${CAMPAIGN}" \
  --source-audit-root "${SOURCE_AUDIT}" \
  --output "${STAGE}/disposition" \
  --expected-plan-artifacts-sha256 "${PLAN_SHA256}" \
  --expected-source-audit-artifacts-sha256 "${SOURCE_AUDIT_SHA256}" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/freeze_external_cohort_development_readiness.py \
  slurm/freeze_gse105127_development_readiness_cpu.sbatch \
  "${VALIDATION}/ARTIFACTS.json" \
  "${PLAN}/ARTIFACTS.json" \
  "${SOURCE_AUDIT}/ARTIFACTS.json" \
  "${CAMPAIGN}/activation/ARTIFACTS.json" \
  "${CAMPAIGN}/terminal/ARTIFACTS.json" \
  "${CAMPAIGN}/cpg_crosswalk/ARTIFACTS.json" \
  "${CAMPAIGN}/rna_consolidated/ARTIFACTS.json" \
  "${STAGE}/disposition/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${VALIDATION_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
validation_sha256 = sys.argv[2]
manifest = verify_frozen_tree(root / "disposition")
receipt = json.loads(
    (root / "disposition/development_readiness.json").read_text(encoding="utf-8")
)
if (
    manifest["metadata"].get("cohort") != "GSE105127"
    or receipt.get("status") != "assay_inputs_ready_development_modeling_blocked"
    or receipt.get("modeling_authorized_by_this_receipt") is not False
    or receipt.get("outcomes_read") is not False
    or receipt.get("fit_or_score_performed") is not False
):
    raise SystemExit("GSE105127 readiness disposition differs")
freeze_tree(
    root,
    {
        "artifact_class": "external_cohort_development_readiness_execution",
        "cohort": "GSE105127",
        "validation_artifacts_sha256": validation_sha256,
        "disposition_artifacts_sha256": sha256_file(
            root / "disposition/ARTIFACTS.json"
        ),
        "modeling_authorized": False,
        "status": receipt["status"],
    },
)
verify_frozen_tree(root)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/disposition/ARTIFACTS.json"
