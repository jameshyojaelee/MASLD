#!/usr/bin/env bash
#SBATCH --job-name=data-check-253
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
PLAN=${ROOT}/executions/model-data-075-21082574/plan
PLAN_SHA256=29566d52311a4204adbfd7ec4501bf07b452ca770ff81f47a3da47fc89b2d717
PRIOR_MODEL_INPUT=${ROOT}/executions/model-data-074-21080832/fixture/model_input
PRIOR_MODEL_INPUT_SHA256=75047dff8fd1e1bc994bf1f915f80c81c87038af5a78c65931fe2faa6c6d2294
CONSOLIDATION=${ROOT}/executions/model-data-079-21083342
# Overridable so a re-run does not require editing this file. The default is
# the original audit execution; pass INDEPENDENT_AUDIT by --export to point
# at a newer one. Matches how this chain already passes INDEPENDENT_VALIDATION_ROOT
# and READINESS_VALIDATION_ROOT.
: "${INDEPENDENT_AUDIT:=${ROOT}/executions/model-data-082-21089146}"
TARGET=${ROOT}/executions/data-check-253-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE268273 readiness disposition.\n' >&2
  exit 1
fi
# The script owns ${STAGE}/disposition: freeze_external_cohort_development_
# readiness.py creates it with mkdir(mode=0o750), no exist_ok, and refuses if
# it already exists. Pre-creating it here made that guard fire on this
# wrapper's own scaffolding -- the same defect that killed the audit step.
install -d -m 0750 "${STAGE}"
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
bash -n slurm/freeze_gse268273_development_readiness_cpu.sbatch
# The original idiom created ${STAGE}/disposition and removed it again here so
# the script could own it. Creating only ${STAGE} above achieves the same
# thing with one step instead of two, so this rmdir is no longer needed.
"${PYTHON}" scripts/freeze_external_cohort_development_readiness.py gse268273 \
  --plan-root "${PLAN}" \
  --prior-model-input-root "${PRIOR_MODEL_INPUT}" \
  --consolidation-root "${CONSOLIDATION}" \
  --independent-audit-root "${INDEPENDENT_AUDIT}" \
  --output "${STAGE}/disposition" \
  --expected-plan-artifacts-sha256 "${PLAN_SHA256}" \
  --expected-prior-model-input-artifacts-sha256 "${PRIOR_MODEL_INPUT_SHA256}" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/freeze_external_cohort_development_readiness.py \
  slurm/freeze_gse268273_development_readiness_cpu.sbatch \
  "${VALIDATION}/ARTIFACTS.json" \
  "${PLAN}/ARTIFACTS.json" \
  "${PRIOR_MODEL_INPUT}/ARTIFACTS.json" \
  "${CONSOLIDATION}/ARTIFACTS.json" \
  "${CONSOLIDATION}/model_input/ARTIFACTS.json" \
  "${INDEPENDENT_AUDIT}/ARTIFACTS.json" \
  "${INDEPENDENT_AUDIT}/audit/ARTIFACTS.json" \
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
    manifest["metadata"].get("cohort") != "GSE268273"
    or receipt.get("status")
    != "raw_source_mechanics_ready_development_training_blocked"
    or receipt.get("modeling_authorized_by_this_receipt") is not False
    or receipt.get("outcomes_read") is not False
    or receipt.get("fit_or_score_performed") is not False
):
    raise SystemExit("GSE268273 readiness disposition differs")
freeze_tree(
    root,
    {
        "artifact_class": "external_cohort_development_readiness_execution",
        "cohort": "GSE268273",
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
