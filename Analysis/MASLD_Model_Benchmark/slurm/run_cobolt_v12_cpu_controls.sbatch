#!/usr/bin/env bash
#SBATCH --job-name=model-training-218
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=18:00:00
#SBATCH --signal=B:USR1@3600
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-218-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-218-%j.err

# Re-fit the 150 CPU control runs under campaign v12.
#
# run_id = canonical_sha256(identity_payload) includes campaign_id and
# code_lock_sha256, so the B6K runtime revision re-keys every run in the campaign,
# including these 150 whose own adapter (rna_atac_classical.py) did not change.
# The completed v9 matrix at cobolt-v9-cpu-baselines--21088479 therefore cannot be
# joined to v12 GPU fits without crossing two source selection records, which the repository
# forbids.  This bundle supersedes it FOR v12 COMPARISON ONLY; the v9 matrix stays
# intact as its own frozen record and is neither deleted nor modified.
#
# Runs on `cpu`, consumes no GPU, and is independent of the B6K GPU bundle, so it
# costs nothing on the critical path.
#
# Required exports: COBOLT_V12_CPU_VALIDATION, COBOLT_V12_CPU_VALIDATION_SHA

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
: "${COBOLT_V12_CPU_VALIDATION:?exact CPU validation artifact is required}"
: "${COBOLT_V12_CPU_VALIDATION_SHA:?exact CPU validation ARTIFACTS SHA-256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v12--3972f43112b8ecc4
CANDIDATE_SHA=9b878eb08b3fc8a6012f473694b778a771da1bdd976f59e2966776f0bda3519f
STABILITY=${ROOT}/executions/cobolt-v12-revision-scope-validation-21111670
STABILITY_SHA=0952f39b35c6d9ea437a32742bfacaf10d3a513e4507d6812b5ca6174e132e0a
EXECUTOR=${ROOT}/scripts/run_cobolt_v12_cpu_baseline_bundle.py
EXECUTOR_SHA=bd710456aff1e307854b68f0edd214e02f181656c151d35851c1d291e1f8e1fa
EXECUTION_ROOT=${ROOT}/candidates/executions/v1-rna-atac-cobolt-5seed-smoke-v12--3972f43112b8ecc4
RECEIPT_ROOT=${ROOT}/executions/cobolt-v12-cpu-control-bundles-shared-root

for digest in "${CANDIDATE_SHA}" "${STABILITY_SHA}" "${EXECUTOR_SHA}" "${COBOLT_V12_CPU_VALIDATION_SHA}"; do
  if [[ ${digest} == *[!0-9a-f]* || ${#digest} -ne 64 ]]; then
    printf 'Invalid SHA-256 binding.\n' >&2
    exit 1
  fi
done
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/model-training-218-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EXECUTOR_SHA}" "${EXECUTOR}" | sha256sum --check --strict
printf '%s  %s\n' "${COBOLT_V12_CPU_VALIDATION_SHA}" "${COBOLT_V12_CPU_VALIDATION}/ARTIFACTS.json" | sha256sum --check --strict
export PYTHONPATH=${ROOT}:${ROOT}/src
python "${EXECUTOR}" \
  --candidate "${CANDIDATE}" \
  --candidate-sha256 "${CANDIDATE_SHA}" \
  --stability "${STABILITY}" \
  --stability-sha256 "${STABILITY_SHA}" \
  --executor-sha256 "${EXECUTOR_SHA}" \
  --workers 12 \
  --allocated-cpus 16 \
  --allocated-memory-gb 200 \
  --execution-root "${EXECUTION_ROOT}" \
  --receipt-root "${RECEIPT_ROOT}" \
  --execution-id "${SLURM_JOB_ID}" \
  --execute
