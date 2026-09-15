#!/usr/bin/env bash
#SBATCH --job-name=model-training-217
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=18:00:00
#SBATCH --signal=B:USR1@3600
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-217-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-217-%j.err

# Re-fit the 150 CPU control runs under campaign v11.
#
# run_id = canonical_sha256(identity_payload) includes campaign_id and
# code_lock_sha256, so the B6K runtime revision re-keys every run in the campaign,
# including these 150 whose own adapter (rna_atac_classical.py) did not change.
# The completed v9 matrix at cobolt-v9-cpu-baselines--21088479 therefore cannot be
# joined to v11 GPU fits without crossing two source selection records, which the repository
# forbids.  This bundle supersedes it FOR v11 COMPARISON ONLY; the v9 matrix stays
# intact as its own frozen record and is neither deleted nor modified.
#
# Runs on `cpu`, consumes no GPU, and is independent of the B6K GPU bundle, so it
# costs nothing on the critical path.
#
# Required exports: COBOLT_V11_CPU_VALIDATION, COBOLT_V11_CPU_VALIDATION_SHA

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
: "${COBOLT_V11_CPU_VALIDATION:?exact CPU validation artifact is required}"
: "${COBOLT_V11_CPU_VALIDATION_SHA:?exact CPU validation ARTIFACTS SHA-256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v11--e17be3dab9ba48cb
CANDIDATE_SHA=dc4f962c841f6daaf50eadcbd476d1ed26eb0a57deada67222fa3ff998c71f47
STABILITY=${ROOT}/executions/cobolt-v11-revision-scope-validation-21109214
STABILITY_SHA=2a3e5c24569c615721519e006c1a61120637dc5697c50153915952478cd6f417
EXECUTOR=${ROOT}/scripts/run_cobolt_v11_cpu_baseline_bundle.py
EXECUTOR_SHA=8567ce270d37db6ab215611e03a695eb5c2fe34afddc357d1dee0c09a4576eeb
EXECUTION_ROOT=${ROOT}/executions/campaign-runs/v1-rna-atac-cobolt-5seed-smoke-v11--e17be3dab9ba48cb--cpu-controls
RECEIPT_ROOT=${ROOT}/executions/cobolt-v11-cpu-control-bundles

for digest in "${CANDIDATE_SHA}" "${STABILITY_SHA}" "${EXECUTOR_SHA}" "${COBOLT_V11_CPU_VALIDATION_SHA}"; do
  if [[ ${digest} == *[!0-9a-f]* || ${#digest} -ne 64 ]]; then
    printf 'Invalid SHA-256 binding.\n' >&2
    exit 1
  fi
done
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/model-training-217-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EXECUTOR_SHA}" "${EXECUTOR}" | sha256sum --check --strict
printf '%s  %s\n' "${COBOLT_V11_CPU_VALIDATION_SHA}" "${COBOLT_V11_CPU_VALIDATION}/ARTIFACTS.json" | sha256sum --check --strict
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
