#!/usr/bin/env bash
#SBATCH --job-name=model-data-050
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PACKAGE_ROOT="${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark"
CONTRACT_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration/candidates/preflight-2026-08-13T224413Z-20206577/contract"
ATLAS="${PROJECT_ROOT}/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
PYTHON_BIN="${PROJECT_ROOT}/.agent/scanpy_env/bin/python"
BUILDER="${PACKAGE_ROOT}/scripts/build_resource_atlas_frozen_screen_50000.py"
UNIT_TEST="tests.unit.test_build_resource_atlas_frozen_screen_50000"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="${PACKAGE_ROOT}/executions/model-data-050-${RUN_STAMP}-${SLURM_JOB_ID}"
FIXTURE="${RUN_DIR}/fixture"
RUN_POINTER="${PACKAGE_ROOT}/executions/model-data-050-${SLURM_JOB_ID}.run_dir"

mkdir "${RUN_DIR}"
if [[ -e "${RUN_POINTER}" ]]; then
  printf 'Refusing to replace run pointer: %s\n' "${RUN_POINTER}" >&2
  exit 1
fi
printf '%s\n' "${RUN_DIR}" > "${RUN_POINTER}"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${PACKAGE_ROOT}/src:${PACKAGE_ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

cd "${PACKAGE_ROOT}"
"${PYTHON_BIN}" -m unittest "${UNIT_TEST}" \
  > "${RUN_DIR}/unit_tests.txt" 2>&1

"${PYTHON_BIN}" "${BUILDER}" \
  --atlas "${ATLAS}" \
  --contract-lock "${CONTRACT_ROOT}/contract_lock.json" \
  --cell-manifest "${CONTRACT_ROOT}/cell_manifest.tsv.gz" \
  --output "${FIXTURE}" \
  --expected-contract-file-sha256 ee67de1c51116285f52cbc2c31a69b468b9ebb585a36800311bdb9cc117ee4ed \
  --expected-contract-lock-sha256 a3730cf2043595dd8318fc20510b4d9ca8b46cf2bd9e5539c6f8424d73bfe12b \
  --expected-manifest-sha256 ef1ffcca70d396aaf9f6f9da0b59bb113ed521aa9f8c989cd12fb194c316805f \
  --expected-raw-counts-sha256 5d405429b69bec62b00ee4a36c19dc124959561acf25d8f5dc384fa8d247c0e0 \
  --expected-cell-order-sha256 f5c8ef5f2747de2649fe3dd8e1ab8a2dcb7ad17223a1e4ee434496904401c15f \
  --expected-gene-order-sha256 a22a5e832dc91a28fc0c01413b835816e543c846c6badcdd03a0a7f8af4f82d7 \
  > "${RUN_DIR}/builder_stdout.txt" 2> "${RUN_DIR}/builder_stderr.txt"

module load python/3.11.5-GCCcore-13.2.0
python - "${FIXTURE}" > "${RUN_DIR}/frozen_tree_verification.json" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

print(json.dumps(verify_frozen_tree(Path(sys.argv[1])), sort_keys=True))
PY
module purge

"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"
sha256sum \
  "${BUILDER}" \
  "${PACKAGE_ROOT}/tests/unit/test_build_resource_atlas_frozen_screen_50000.py" \
  "${PACKAGE_ROOT}/slurm/build_resource_atlas_frozen_screen_50000.sbatch" \
  > "${RUN_DIR}/source_sha256.txt"
sha256sum \
  "${CONTRACT_ROOT}/contract_lock.json" \
  "${CONTRACT_ROOT}/cell_manifest.tsv.gz" \
  > "${RUN_DIR}/input_sha256.txt"
sha256sum "${FIXTURE}/ARTIFACTS.json" "${FIXTURE}/COMPLETE"
