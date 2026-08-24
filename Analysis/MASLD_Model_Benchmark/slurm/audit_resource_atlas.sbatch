#!/usr/bin/env bash
#SBATCH --job-name=masld-atlas-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PACKAGE_ROOT="${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="${PACKAGE_ROOT}/executions/resource-atlas-activation-audit-${RUN_STAMP}-${SLURM_JOB_ID}"

mkdir -p "${RUN_DIR}"
printf '%s\n' "${RUN_DIR}" > "${PACKAGE_ROOT}/executions/masld-atlas-audit-${SLURM_JOB_ID}.run_dir"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

"${PYTHON_BIN}" "${PACKAGE_ROOT}/src/masld_bench/activation_audit.py" \
  --atlas "${PROJECT_ROOT}/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad" \
  --contract-lock "${PROJECT_ROOT}/Analysis/SingleCell/continual_integration/candidates/preflight-2026-08-13T224413Z-20206577/contract/contract_lock.json" \
  --output "${RUN_DIR}/audit"

"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"
sha256sum \
  "${PACKAGE_ROOT}/src/masld_bench/activation_audit.py" \
  "${PACKAGE_ROOT}/slurm/audit_resource_atlas.sbatch" \
  > "${RUN_DIR}/source_sha256.txt"
