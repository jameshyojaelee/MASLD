#!/usr/bin/env bash
#SBATCH --job-name=model-data-082
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${INDEPENDENT_VALIDATION_ROOT:?INDEPENDENT_VALIDATION_ROOT must name the frozen independent validation}"
: "${INDEPENDENT_VALIDATION_SHA256:?INDEPENDENT_VALIDATION_SHA256 must pin that ARTIFACTS.json}"

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
CONSOLIDATION=${ROOT}/executions/model-data-079-21083342
QUANT_STATE=${ROOT}/executions/resumable-model-data-078/fb28496669135258d3a6e9895299cb2757453056038aa2b4fd6fba05516f8d3a
PRIOR_MODEL_INPUT=${ROOT}/executions/model-data-074-21080832/fixture/model_input
PRIOR_MODEL_INPUT_SHA256=75047dff8fd1e1bc994bf1f915f80c81c87038af5a78c65931fe2faa6c6d2294
PLAN=${ROOT}/executions/model-data-075-21082574/plan
PLAN_SHA256=29566d52311a4204adbfd7ec4501bf07b452ca770ff81f47a3da47fc89b2d717
REFERENCE=${ROOT}/executions/model-data-077-21083340/artifact
REFERENCE_SHA256=8c4282a6070f6eb99a8158ac4f496a4715ee47042d243e13b69ce4b0f1915672
ADMISSION_EXECUTION=${ROOT}/executions/model-data-081-21088737
ADMISSION_SOURCE_SHA256=db6bc3538749a980bd8606b131b906fb30d2e82d7498df057c4986db9eb72785
ADMISSION_WRAPPER_SHA256=4dd743bfa85df4fb8b1c690255623f7fdbbe5c845955198e0190d2ff19a8f6d5
CORE_VALIDATION_WRAPPER_SHA256=ba7b32e091658b11f213e0a65219965030814dc7147d46f41debef1d29c6ebbb
PRODUCTION_VALIDATION_WRAPPER_SHA256=c2bc3ef6fd5c8c283c7f604aedd889417b2008a04ab2e7f1061572d41a11c0f2
CORE_VALIDATION_ARTIFACTS_SHA256=1b28d012e5f593a7ccfdd073d2872deb6f47d7b247ffe8b5db4026910fb110cd
WRAPPER_VALIDATION_ARTIFACTS_SHA256=865ea6a736d6c004c2db1a7ba8ff81d40e3a73a69f0cad51375f319b74661924
TARGET=${ROOT}/executions/model-data-082-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite independent GSE268273 source audit.\n' >&2
  exit 1
fi
# Create the stage only. The audit script owns ${STAGE}/audit and creates it
# itself at scripts/audit_gse268273_raw_source_integrity_unexposed.py:476 with
# output.mkdir(mode=0o750) and no exist_ok, then refuses at :254 if it already
# exists. Pre-creating it here claimed ownership a second time and made that
# guard fire on this wrapper's own scaffolding, so the job could never succeed
# on any run. Do not add /audit back.
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${INDEPENDENT_VALIDATION_SHA256}" "${INDEPENDENT_VALIDATION_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PRIOR_MODEL_INPUT_SHA256}" "${PRIOR_MODEL_INPUT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PLAN_SHA256}" "${PLAN}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${REFERENCE_SHA256}" "${REFERENCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${ADMISSION_SOURCE_SHA256}" "${ROOT}/scripts/admit_gse268273_raw_rsem_source_integrity.py" \
  | sha256sum --check --strict
printf '%s  %s\n' "${ADMISSION_WRAPPER_SHA256}" "${ROOT}/slurm/admit_gse268273_raw_source_integrity_cpu.sbatch" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CORE_VALIDATION_WRAPPER_SHA256}" "${ROOT}/slurm/validate_gse268273_raw_source_integrity_admission_cpu.sbatch" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PRODUCTION_VALIDATION_WRAPPER_SHA256}" "${ROOT}/slurm/validate_gse268273_raw_source_integrity_production_cpu.sbatch" \
  | sha256sum --check --strict
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
bash -n slurm/audit_gse268273_raw_source_integrity_unexposed_cpu.sbatch
{
  printf 'JobID\tJobName\tPartition\tQOS\tState\tExitCode\tSubmitLine\n'
  sacct -X -j 21088737 -n -P \
    -o JobID,JobName,Partition,QOS,State,ExitCode,SubmitLine%1000 | tr '|' '\t'
} > "${STAGE}/admission_submission.tsv"
"${PYTHON}" scripts/audit_gse268273_raw_source_integrity_unexposed.py \
  --consolidation-root "${CONSOLIDATION}" \
  --quant-state-root "${QUANT_STATE}" \
  --prior-model-input-root "${PRIOR_MODEL_INPUT}" \
  --plan-root "${PLAN}" \
  --reference-root "${REFERENCE}" \
  --admission-execution-root "${ADMISSION_EXECUTION}" \
  --admission-source "${ROOT}/scripts/admit_gse268273_raw_rsem_source_integrity.py" \
  --admission-wrapper "${ROOT}/slurm/admit_gse268273_raw_source_integrity_cpu.sbatch" \
  --core-validation-wrapper "${ROOT}/slurm/validate_gse268273_raw_source_integrity_admission_cpu.sbatch" \
  --production-validation-wrapper "${ROOT}/slurm/validate_gse268273_raw_source_integrity_production_cpu.sbatch" \
  --submission-record "${STAGE}/admission_submission.tsv" \
  --output "${STAGE}/audit" \
  --expected-prior-model-input-artifacts-sha256 "${PRIOR_MODEL_INPUT_SHA256}" \
  --expected-plan-artifacts-sha256 "${PLAN_SHA256}" \
  --expected-reference-artifacts-sha256 "${REFERENCE_SHA256}" \
  --expected-admission-source-sha256 "${ADMISSION_SOURCE_SHA256}" \
  --expected-admission-wrapper-sha256 "${ADMISSION_WRAPPER_SHA256}" \
  --expected-core-validation-wrapper-sha256 "${CORE_VALIDATION_WRAPPER_SHA256}" \
  --expected-production-validation-wrapper-sha256 "${PRODUCTION_VALIDATION_WRAPPER_SHA256}" \
  --expected-core-validation-artifacts-sha256 "${CORE_VALIDATION_ARTIFACTS_SHA256}" \
  --expected-wrapper-validation-artifacts-sha256 "${WRAPPER_VALIDATION_ARTIFACTS_SHA256}" \
  --expected-admission-job-id 21088737 \
  --expected-participants 109 \
  --expected-genes 14078 \
  --expected-technical-runs 824 \
  --expected-fastq-files 968 \
  --expected-fastq-bytes 512881099727 \
  --expected-bundles 8 \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/audit_gse268273_raw_source_integrity_unexposed.py \
  slurm/audit_gse268273_raw_source_integrity_unexposed_cpu.sbatch \
  scripts/admit_gse268273_raw_rsem_source_integrity.py \
  slurm/admit_gse268273_raw_source_integrity_cpu.sbatch \
  slurm/validate_gse268273_raw_source_integrity_admission_cpu.sbatch \
  slurm/validate_gse268273_raw_source_integrity_production_cpu.sbatch \
  "${CONSOLIDATION}/ARTIFACTS.json" \
  "${CONSOLIDATION}/model_input/ARTIFACTS.json" \
  "${PRIOR_MODEL_INPUT}/ARTIFACTS.json" \
  "${PLAN}/ARTIFACTS.json" \
  "${REFERENCE}/ARTIFACTS.json" \
  "${ADMISSION_EXECUTION}/ARTIFACTS.json" \
  "${INDEPENDENT_VALIDATION_ROOT}/ARTIFACTS.json" \
  "${STAGE}/admission_submission.tsv" \
  "${STAGE}/audit/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" "${INDEPENDENT_VALIDATION_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
validation_sha256 = sys.argv[2]
manifest = verify_frozen_tree(root / "audit")
receipt = json.loads(
    (root / "audit/independent_source_integrity_audit.json").read_text(encoding="utf-8")
)
if (
    manifest["metadata"].get("eligible_scope") != "source_integrity_mechanics_only"
    or receipt.get("source_integrity_mechanics_eligible") is not True
    or receipt.get("modeling_or_scoring_authorized_by_this_receipt") is not False
    or receipt.get("sealed_evaluation_authorized_by_this_receipt") is not False
    or receipt.get("clean_or_sealed_champion_eligible_by_this_receipt") is not False
):
    raise SystemExit("independent GSE268273 source-only audit firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse268273_independent_raw_source_integrity_audit_execution",
        "validation_artifacts_sha256": validation_sha256,
        "source_integrity_mechanics_eligible": True,
        "eligible_scope": "source_integrity_mechanics_only",
        "modeling_or_scoring_authorized": False,
        "clean_or_sealed_champion_eligible": False,
        "status": "passed_source_integrity_mechanics_only",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/audit/ARTIFACTS.json"
