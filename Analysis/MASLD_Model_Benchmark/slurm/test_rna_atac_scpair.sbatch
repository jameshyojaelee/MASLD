#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-scpair-test
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
PYTHON_BIN="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python"
ENV_LOCK="${ROOT}/executions/environment-locks/gpu_rna_atac_torch_smoke--c2308fc3b7524081/environment.lock.json"
DATA="${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset/gse296875_rna_atac_smoke_1000.h5"
OUTPUT="${ROOT}/executions/tests/rna-atac-scpair-e2e-${SLURM_JOB_ID}"

export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${ROOT}/src"
export CUBLAS_WORKSPACE_CONFIG=:4096:8
cd "${ROOT}"

"${PYTHON_BIN}" -m py_compile \
  src/masld_bench/adapters/rna_atac_scpair.py \
  tests/scientific/test_rna_atac_scpair_e2e.py \
  tests/unit/test_rna_atac_scpair.py
"${PYTHON_BIN}" -m unittest tests.unit.test_rna_atac_scpair -v
"${PYTHON_BIN}" tests/scientific/test_rna_atac_scpair_e2e.py \
  --environment-lock "${ENV_LOCK}" \
  --data "${DATA}" \
  --output "${OUTPUT}"

"${PYTHON_BIN}" -m pip freeze > "${OUTPUT}/pip_freeze.txt"
sha256sum \
  "${ROOT}/src/masld_bench/adapters/rna_atac_scpair.py" \
  "${ROOT}/tests/scientific/test_rna_atac_scpair_e2e.py" \
  "${ROOT}/tests/unit/test_rna_atac_scpair.py" \
  "${ROOT}/slurm/test_rna_atac_scpair.sbatch" \
  > "${OUTPUT}/source_sha256.txt"
