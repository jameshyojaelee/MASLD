#!/bin/bash
#SBATCH --job-name=ffreeze
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/resource_f_five/ffreeze_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/resource_f_five/ffreeze_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/manuscript/resource_f_five"
CANDIDATE_ROOT="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10"
BG_ROOT="${PROJECT_ROOT}/results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T195243Z"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
RUNTIME_VERIFIER="${BG_ROOT}/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/analysis_runtime_contract.py"
RUNTIME_SAFE_IO="${BG_ROOT}/source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/safe_io.py"
RNASEQ_PYTHON="${RNASEQ_ENV}/bin/python3.10"
EXPECTED_PYTHON_SHA256="2e76c2d5e89880deb45cf75abfb81a10029f1d9cbbb155686e036720f72f7502"
EXPECTED_RUNTIME_VERIFIER_SHA256="7921e8468f6bc08cd508c996224eef0f3c7e07f7ba09c59dd78b77ce2adc7342"
EXPECTED_RUNTIME_SAFE_IO_SHA256="2b98c1b94d4392c96055ef9f92fab020bb314fc17fa58955fe13147212d72cc5"
REVIEWED_MANIFEST="${SCRIPT_DIR}/reviewed_code.sha256"

cd "${PROJECT_ROOT}"

for pinned_file in "${RNASEQ_PYTHON}" "${RUNTIME_VERIFIER}" "${RUNTIME_SAFE_IO}"; do
  if [[ ! -f "${pinned_file}" || -L "${pinned_file}" ]]; then
    echo "Pinned runtime file is missing or symlinked: ${pinned_file}" >&2
    exit 1
  fi
done
if [[ "$(readlink -f -- "${RNASEQ_PYTHON}")" != "${RNASEQ_PYTHON}" ]]; then
  echo "Pinned Python resolved-path drift" >&2
  exit 1
fi
if [[ "$(sha256sum -- "${RNASEQ_PYTHON}" | cut -d " " -f1)" != "${EXPECTED_PYTHON_SHA256}" ]]; then
  echo "Pinned Python SHA256 drift" >&2
  exit 1
fi
if [[ "$(sha256sum -- "${RUNTIME_VERIFIER}" | cut -d " " -f1)" != "${EXPECTED_RUNTIME_VERIFIER_SHA256}" ]]; then
  echo "Pinned runtime verifier SHA256 drift" >&2
  exit 1
fi
if [[ "$(sha256sum -- "${RUNTIME_SAFE_IO}" | cut -d " " -f1)" != "${EXPECTED_RUNTIME_SAFE_IO_SHA256}" ]]; then
  echo "Pinned runtime safe-I/O SHA256 drift" >&2
  exit 1
fi

unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
PYCACHE_PARENT="${SLURM_TMPDIR:-/tmp}"
if [[ ! -d "${PYCACHE_PARENT}" || -L "${PYCACHE_PARENT}" ]]; then
  echo "Invalid private Python-cache parent: ${PYCACHE_PARENT}" >&2
  exit 1
fi
RUNTIME_PYCACHE=$(mktemp -d "${PYCACHE_PARENT}/resource-f-five-pycache.XXXXXX")
if [[ ! -d "${RUNTIME_PYCACHE}" || -L "${RUNTIME_PYCACHE}" ]]; then
  echo "Failed to create private empty Python-cache root" >&2
  exit 1
fi
export PYTHONPYCACHEPREFIX="${RUNTIME_PYCACHE}"
export PYTHONDONTWRITEBYTECODE=1

if [[ -e "${CANDIDATE_ROOT}" || -L "${CANDIDATE_ROOT}" ]]; then
  echo "Refusing existing Resource candidate: ${CANDIDATE_ROOT}" >&2
  exit 1
fi

: "${MASLD_REVIEWED_CODE_MANIFEST_SHA256:?submit with the independently reviewed manifest SHA256}"
if [[ ! "${MASLD_REVIEWED_CODE_MANIFEST_SHA256}" =~ ^[0-9a-f]{64}$ ]]; then
  echo "Invalid reviewed-code manifest SHA256" >&2
  exit 1
fi
if [[ "$(sha256sum "${REVIEWED_MANIFEST}" | awk '{print $1}')" != "${MASLD_REVIEWED_CODE_MANIFEST_SHA256}" ]]; then
  echo "Reviewed-code manifest identity drift" >&2
  exit 1
fi
sha256sum --strict -c "${REVIEWED_MANIFEST}"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONDONTWRITEBYTECODE=1
export CONDA_PREFIX="${RNASEQ_ENV}"
export PATH="${RNASEQ_ENV}/bin:/gpfs/commons/home/jameslee/.local/bin:${PATH}"
verify_frozen_runtime() (
  # The accepted runtime was sealed with user-site package visibility. Keep
  # Resource Python isolated, but let micromamba observe that sealed surface.
  unset PYTHONNOUSERSITE
  "${RNASEQ_PYTHON}" -S "${RUNTIME_VERIFIER}" verify --run-root "$1"
)

verify_frozen_runtime "${BG_ROOT}"
"${RNASEQ_PYTHON}" -S "${SCRIPT_DIR}/prepare_candidate.py"

test -f "${CANDIDATE_ROOT}/BASE_SNAPSHOT_COMPLETE.json"
test -f "${CANDIDATE_ROOT}/manifests/source_selection.tsv"
test -f "${CANDIDATE_ROOT}/inputs/BG001-DECISION/arms/F_five/ARM_COMPLETE.json"
test -f "${CANDIDATE_ROOT}/inputs/BG001-DECISION/arms/F_legacy/ARM_COMPLETE.json"
"${RNASEQ_PYTHON}" -S "${CANDIDATE_ROOT}/code/verify_candidate_snapshot.py"
echo "BASE_SNAPSHOT_COMPLETE: ${CANDIDATE_ROOT}"
