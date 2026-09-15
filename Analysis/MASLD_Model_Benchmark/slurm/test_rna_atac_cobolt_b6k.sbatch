#!/usr/bin/env bash
#SBATCH --job-name=model-probe-013
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:b6k:1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Cobolt adapter test suite under the B6K cu128 runtime, ahead of the v11 freeze.
# Mirrors slurm/test_rna_atac_cobolt.sbatch but binds the new runtime and the new
# environment pin, and adds the allowlist fail-closed tests.
#
# Submitted directly rather than through the GPU dependency dispatcher, which
# hardcodes gpu:l40s:1 and cannot route B6K.  Job ids from this wrapper never
# enter the dispatcher roster, so job 21100574 is unaffected.
#
# Produces no campaign fit and computes no benchmark metric.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME_TREE="${ROOT}/executions/environments/cobolt-native-torch2.8.0-cu128-21107432"
RUNTIME_TREE_SHA256=e7e56c0ee207843b78dd9bbb1e3619fc965fc7836a4d3b0131c58b6f4e58049c
PYTHON="${RUNTIME_TREE}/env/bin/python"
ENVIRONMENT="${ROOT}/executions/environment-locks/gpu_rna_atac_torch_b6k--7b1948ff69fb8a46/environment.lock.json"
ENVIRONMENT_SHA256=7b1948ff69fb8a46831371424df0f02c897767c7a061af832d513f2d892f50ca
DATA="${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset/gse296875_rna_atac_smoke_1000.h5"
UPSTREAM="${ROOT}/executions/upstream-acquisitions/cobolt-v1.0.1-21060709/source/epurdom-cobolt-cf5a448/cobolt/model/coboltmodel.py"
UPSTREAM_SHA256=a4f5c3a45ca4069e12ffa05d3efbfcd8a0a636da7515d7d4e9174a914f3f2311
OUT="${ROOT}/executions/cobolt-b6k-adapter-test-${SLURM_JOB_ID}"

cd "${ROOT}"
sha256sum --check --strict <<EOF
${RUNTIME_TREE_SHA256}  ${RUNTIME_TREE}/ARTIFACTS.json
${ENVIRONMENT_SHA256}  ${ENVIRONMENT}
${UPSTREAM_SHA256}  ${UPSTREAM}
EOF

nvidia-smi --query-gpu=name --format=csv,noheader

# Tests that load the adapter standalone run under the B6K runtime itself.
PYTHONPATH=src "${PYTHON}" -m unittest -v \
  tests.unit.test_rna_atac_cobolt \
  tests.unit.test_rna_atac_regularized

# test_rna_atac_classical and test_rna_atac_cobolt_evaluator import the
# masld_bench package, whose contracts.py needs stdlib tomllib (python 3.11+).
# Both GPU runtimes are python 3.10, so these are control-plane tests and run
# under the module python, exactly as the rest of the repository runs them.
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_rna_atac_classical \
  tests.unit.test_rna_atac_cobolt_evaluator

PYTHONPATH=src "${PYTHON}" tests/scientific/test_rna_atac_cobolt_e2e.py \
  --environment-lock "${ENVIRONMENT}" \
  --data "${DATA}" \
  --upstream-source "${UPSTREAM}" \
  --output "${OUT}"

printf 'Cobolt B6K adapter test suite passed: %s\n' "${OUT}"
