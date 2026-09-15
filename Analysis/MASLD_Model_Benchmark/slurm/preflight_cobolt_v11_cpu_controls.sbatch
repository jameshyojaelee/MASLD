#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-354
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-354-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-check-354-%j.err

# Preflight for the v11 CPU control re-run. Verifies every shared input and the
# 150-run selection WITHOUT executing any run, and freezes the validation file
# that the execute wrapper then binds by SHA-256.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v11--e17be3dab9ba48cb
CANDIDATE_SHA=dc4f962c841f6daaf50eadcbd476d1ed26eb0a57deada67222fa3ff998c71f47
STABILITY=${ROOT}/executions/cobolt-v11-revision-scope-validation-21109214
STABILITY_SHA=2a3e5c24569c615721519e006c1a61120637dc5697c50153915952478cd6f417
EXECUTOR=${ROOT}/scripts/run_cobolt_v11_cpu_baseline_bundle.py
EXECUTOR_SHA=8567ce270d37db6ab215611e03a695eb5c2fe34afddc357d1dee0c09a4576eeb
RECEIPT_ROOT=${ROOT}/executions/cobolt-v11-cpu-control-preflights

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/cobolt-v11-preflight-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EXECUTOR_SHA}" "${EXECUTOR}" | sha256sum --check --strict
export PYTHONPATH=${ROOT}:${ROOT}/src
python "${EXECUTOR}" \
  --candidate "${CANDIDATE}" \
  --candidate-sha256 "${CANDIDATE_SHA}" \
  --stability "${STABILITY}" \
  --stability-sha256 "${STABILITY_SHA}" \
  --executor-sha256 "${EXECUTOR_SHA}" \
  --allocated-cpus 2 \
  --allocated-memory-gb 32 \
  --workers 1 \
  --receipt-root "${RECEIPT_ROOT}" \
  --execution-id "${SLURM_JOB_ID}" \
  --preflight
