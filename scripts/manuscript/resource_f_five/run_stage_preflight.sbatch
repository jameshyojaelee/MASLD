#!/bin/bash
#SBATCH --job-name=fstagepre
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/resource_f_five/fstagepre_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/resource_f_five/fstagepre_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE_ROOT="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10"
CODE_ROOT="${CANDIDATE_ROOT}/code"
PREFLIGHT_ROOT="${CANDIDATE_ROOT}/workstreams/BULK-STAGE/PREFLIGHT-v1"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
RUNTIME_ROOT="${CANDIDATE_ROOT}/inputs/BG001-DECISION"
RUNTIME_VERIFIER="${CODE_ROOT}/frozen_bg001_runtime/analysis_runtime_contract.py"
RUNTIME_SAFE_IO="${CODE_ROOT}/frozen_bg001_runtime/safe_io.py"
RNASEQ_PYTHON="${RNASEQ_ENV}/bin/python3.10"
EXPECTED_PYTHON_SHA256="2e76c2d5e89880deb45cf75abfb81a10029f1d9cbbb155686e036720f72f7502"
EXPECTED_RUNTIME_VERIFIER_SHA256="7921e8468f6bc08cd508c996224eef0f3c7e07f7ba09c59dd78b77ce2adc7342"
EXPECTED_RUNTIME_SAFE_IO_SHA256="2b98c1b94d4392c96055ef9f92fab020bb314fc17fa58955fe13147212d72cc5"

test -f "${CANDIDATE_ROOT}/BASE_SNAPSHOT_COMPLETE.json"

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
test -f "${CODE_ROOT}/preflight_stage_substrate.R"
test -f "${CODE_ROOT}/validate_stage_preflight.py"
test -f "${CODE_ROOT}/validate_stage_rds_membership.R"

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export MASLD_RESOURCE_CANDIDATE_ROOT="${CANDIDATE_ROOT}"
export MASLD_RSCRIPT_PATH="${RNASEQ_ENV}/bin/Rscript"
export MASLD_PYTHON_PATH="${RNASEQ_PYTHON}"
export PYTHONDONTWRITEBYTECODE=1
export CONDA_PREFIX="${RNASEQ_ENV}"
export PATH="${RNASEQ_ENV}/bin:/gpfs/commons/home/jameslee/.local/bin:${PATH}"
export LANG="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
verify_frozen_runtime() (
  # The accepted runtime was sealed with user-site package visibility. Keep
  # Resource Python isolated, but let micromamba observe that sealed surface.
  unset PYTHONNOUSERSITE
  "${RNASEQ_PYTHON}" -S "${RUNTIME_VERIFIER}" verify --run-root "$1"
)


"${RNASEQ_PYTHON}" -S "${CODE_ROOT}/verify_candidate_snapshot.py"
verify_frozen_runtime "${RUNTIME_ROOT}"
if [[ ! -e "${PREFLIGHT_ROOT}" && ! -L "${PREFLIGHT_ROOT}" ]]; then
  "${RNASEQ_ENV}/bin/Rscript" --vanilla "${CODE_ROOT}/preflight_stage_substrate.R"
elif [[ -d "${PREFLIGHT_ROOT}" && ! -L "${PREFLIGHT_ROOT}" ]]; then
  echo "Existing stage preflight found; entering validation recovery path"
else
  echo "Invalid stage preflight path: ${PREFLIGHT_ROOT}" >&2
  exit 1
fi
"${RNASEQ_PYTHON}" -S "${CODE_ROOT}/validate_stage_preflight.py"
verify_frozen_runtime "${RUNTIME_ROOT}"
"${RNASEQ_PYTHON}" -S "${CODE_ROOT}/verify_candidate_snapshot.py"
"${RNASEQ_PYTHON}" -S "${CODE_ROOT}/validate_stage_preflight.py" --verify-only

test -f "${PREFLIGHT_ROOT}/STAGE_PREFLIGHT_READY.json"
echo "STAGE_PREFLIGHT_READY: ${PREFLIGHT_ROOT}/STAGE_PREFLIGHT_READY.json"
