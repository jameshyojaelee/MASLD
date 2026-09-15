#!/usr/bin/env bash
#SBATCH --job-name=model-validate-209
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
CORE_VALIDATION=${ROOT}/executions/model-validate-208-21088644
CORE_VALIDATION_SHA256=1b28d012e5f593a7ccfdd073d2872deb6f47d7b247ffe8b5db4026910fb110cd
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/model-validate-209-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE268273 production admission validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${CORE_VALIDATION_SHA256}" "${CORE_VALIDATION}/ARTIFACTS.json" \
  | sha256sum --check --strict
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
bash -n slurm/admit_gse268273_raw_source_integrity_cpu.sbatch
"${PYTHON}" - <<'PY' > "${STAGE}/validation/contract.txt"
from pathlib import Path

path = Path("slurm/admit_gse268273_raw_source_integrity_cpu.sbatch")
source = path.read_text(encoding="utf-8")
required = (
    "#SBATCH --partition=cpu",
    "#SBATCH --qos=nslab",
    "#SBATCH --cpus-per-task=2",
    "#SBATCH --mem=16G",
    "#SBATCH --time=02:00:00",
    "CONSOLIDATION=${ROOT}/executions/model-data-079-21083342",
    "PRIOR_MODEL_INPUT_SHA256=75047dff8fd1e1bc994bf1f915f80c81c87038af5a78c65931fe2faa6c6d2294",
    "PLAN_SHA256=29566d52311a4204adbfd7ec4501bf07b452ca770ff81f47a3da47fc89b2d717",
    "REFERENCE_SHA256=8c4282a6070f6eb99a8158ac4f496a4715ee47042d243e13b69ce4b0f1915672",
    "INCIDENT_SHA256=c104ee308d3762cb035134cdb790dbb752e5a8873d1196566420179449fe06df",
    "CORE_VALIDATION_SHA256=1b28d012e5f593a7ccfdd073d2872deb6f47d7b247ffe8b5db4026910fb110cd",
    "WRAPPER_VALIDATION_ROOT must name the frozen wrapper validation",
    "--expected-participants 109",
    "--expected-genes 14078",
    "--expected-technical-runs 824",
    "--expected-fastq-bytes 512881099727",
    '"source_integrity_only": True',
    '"clean_or_sealed_champion_eligible": False',
)
for value in required:
    if value not in source:
        raise SystemExit(f"missing production admission contract: {value}")
for value in (
    "#SBATCH --array",
    "#SBATCH --partition=gpu",
    "evaluator_only",
    '"outcomes.tsv"',
    "gse268273_fibrosis_ood_transfer_task",
    " sbatch ",
):
    if value in source:
        raise SystemExit(f"forbidden production admission surface: {value}")
print("production source-integrity admission contract passed")
PY
sha256sum \
  scripts/admit_gse268273_raw_rsem_source_integrity.py \
  tests/unit/test_admit_gse268273_raw_rsem_source_integrity.py \
  slurm/admit_gse268273_raw_source_integrity_cpu.sbatch \
  slurm/validate_gse268273_raw_source_integrity_production_cpu.sbatch \
  "${CORE_VALIDATION}/ARTIFACTS.json" \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "gse268273_raw_source_integrity_production_validation",
        "core_validation_artifacts_sha256": (
            "1b28d012e5f593a7ccfdd073d2872deb6f47d7b247ffe8b5db4026910fb110cd"
        ),
        "source_integrity_only": True,
        "clean_or_sealed_champion_eligible": False,
        "independent_unexposed_rederivation_or_audit_required": True,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
