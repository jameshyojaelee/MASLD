#!/usr/bin/env bash
#SBATCH --job-name=model-data-081
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${WRAPPER_VALIDATION_ROOT:?WRAPPER_VALIDATION_ROOT must name the frozen wrapper validation}"
: "${WRAPPER_VALIDATION_SHA256:?WRAPPER_VALIDATION_SHA256 must pin that ARTIFACTS.json}"

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CONSOLIDATION=${ROOT}/executions/model-data-079-21083342
PRIOR_MODEL_INPUT=${ROOT}/executions/model-data-074-21080832/fixture/model_input
PRIOR_MODEL_INPUT_SHA256=75047dff8fd1e1bc994bf1f915f80c81c87038af5a78c65931fe2faa6c6d2294
PLAN=${ROOT}/executions/model-data-075-21082574/plan
PLAN_SHA256=29566d52311a4204adbfd7ec4501bf07b452ca770ff81f47a3da47fc89b2d717
REFERENCE=${ROOT}/executions/model-data-077-21083340/artifact
REFERENCE_SHA256=8c4282a6070f6eb99a8158ac4f496a4715ee47042d243e13b69ce4b0f1915672
INCIDENT=${ROOT}/executions/gse268273-evaluator-receipt-incident-20260824
INCIDENT_SHA256=c104ee308d3762cb035134cdb790dbb752e5a8873d1196566420179449fe06df
CORE_VALIDATION=${ROOT}/executions/model-validate-208-21088644
CORE_VALIDATION_SHA256=1b28d012e5f593a7ccfdd073d2872deb6f47d7b247ffe8b5db4026910fb110cd
TARGET=${ROOT}/executions/model-data-081-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE268273 source admission artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PRIOR_MODEL_INPUT_SHA256}" "${PRIOR_MODEL_INPUT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PLAN_SHA256}" "${PLAN}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${REFERENCE_SHA256}" "${REFERENCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${INCIDENT_SHA256}" "${INCIDENT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CORE_VALIDATION_SHA256}" "${CORE_VALIDATION}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${WRAPPER_VALIDATION_SHA256}" "${WRAPPER_VALIDATION_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
bash -n slurm/admit_gse268273_raw_source_integrity_cpu.sbatch
"${PYTHON}" scripts/admit_gse268273_raw_rsem_source_integrity.py \
  --consolidation-root "${CONSOLIDATION}" \
  --prior-model-input-root "${PRIOR_MODEL_INPUT}" \
  --plan-root "${PLAN}" \
  --reference-root "${REFERENCE}" \
  --incident-root "${INCIDENT}" \
  --output "${STAGE}/admission" \
  --expected-prior-model-input-artifacts-sha256 "${PRIOR_MODEL_INPUT_SHA256}" \
  --expected-plan-artifacts-sha256 "${PLAN_SHA256}" \
  --expected-reference-artifacts-sha256 "${REFERENCE_SHA256}" \
  --expected-incident-artifacts-sha256 "${INCIDENT_SHA256}" \
  --expected-participants 109 \
  --expected-genes 14078 \
  --expected-technical-runs 824 \
  --expected-fastq-bytes 512881099727 \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/admit_gse268273_raw_rsem_source_integrity.py \
  slurm/admit_gse268273_raw_source_integrity_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${CONSOLIDATION}/ARTIFACTS.json" \
  "${CONSOLIDATION}/model_input/ARTIFACTS.json" \
  "${PRIOR_MODEL_INPUT}/ARTIFACTS.json" \
  "${PLAN}/ARTIFACTS.json" \
  "${REFERENCE}/ARTIFACTS.json" \
  "${INCIDENT}/ARTIFACTS.json" \
  "${CORE_VALIDATION}/ARTIFACTS.json" \
  "${WRAPPER_VALIDATION_ROOT}/ARTIFACTS.json" \
  "${STAGE}/admission/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - \
  "${STAGE}" \
  "${CORE_VALIDATION_SHA256}" \
  "${WRAPPER_VALIDATION_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
core_validation_sha256 = sys.argv[2]
wrapper_validation_sha256 = sys.argv[3]
manifest = verify_frozen_tree(root / "admission")
receipt = json.loads(
    (root / "admission/source_integrity_admission.json").read_text(encoding="utf-8")
)
if (
    manifest["metadata"].get("clean_or_sealed_champion_eligible") is not False
    or receipt.get("source_integrity_only") is not True
    or receipt.get("clean_or_sealed_champion_eligible") is not False
    or receipt.get("independent_unexposed_rederivation_or_audit_required") is not True
    or receipt.get("labels_read") is not False
    or receipt.get("clinical_covariates_read") is not False
    or receipt.get("evaluation_data_opened") is not False
    or receipt.get("fit_or_score_performed") is not False
):
    raise SystemExit("GSE268273 source-integrity-only firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse268273_raw_source_integrity_admission_execution",
        "core_validation_artifacts_sha256": core_validation_sha256,
        "wrapper_validation_artifacts_sha256": wrapper_validation_sha256,
        "source_integrity_only": True,
        "clean_or_sealed_champion_eligible": False,
        "independent_unexposed_rederivation_or_audit_required": True,
        "status": "passed_unexposed_audit_required",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/admission/ARTIFACTS.json"
