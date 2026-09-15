#!/usr/bin/env bash
#SBATCH --job-name=model-training-213
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=18:00:00
#SBATCH --signal=B:USR1@3600
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-213-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-213-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${COBOLT_V9_CANDIDATE:?set COBOLT_V9_CANDIDATE to the exact frozen production candidate}
CANDIDATE_SHA=${COBOLT_V9_CANDIDATE_SHA:?set COBOLT_V9_CANDIDATE_SHA to its ARTIFACTS.json sha256}
STABILITY=${COBOLT_V9_STABILITY:?set COBOLT_V9_STABILITY to the exact frozen stability artifact}
STABILITY_SHA=${COBOLT_V9_STABILITY_SHA:?set COBOLT_V9_STABILITY_SHA to its ARTIFACTS.json sha256}
EXECUTOR=${ROOT}/scripts/run_cobolt_v9_cpu_baseline_bundle.py
EXECUTOR_SHA=${COBOLT_V9_CPU_EXECUTOR_SHA:?set COBOLT_V9_CPU_EXECUTOR_SHA to the reviewed executor sha256}
EXECUTION_ROOT=${ROOT}/executions/campaign-runs/${CANDIDATE##*/}--cpu-baselines
RECEIPT_ROOT=${ROOT}/executions/cobolt-v9-cpu-baseline-bundles

for digest in "${CANDIDATE_SHA}" "${STABILITY_SHA}" "${EXECUTOR_SHA}"; do
  if [[ "${digest}" == *[!0-9a-f]* || ${#digest} -ne 64 ]]; then
    printf 'Invalid SHA-256 binding\n' >&2
    exit 1
  fi
done
if [[ "${CANDIDATE}" != "${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--7031e0dfbf76bc4d" ]]; then
  printf 'Invalid v9 candidate path\n' >&2
  exit 1
fi

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/model-training-213-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EXECUTOR_SHA}" "${EXECUTOR}" | sha256sum --check --strict

PYTHONPATH=src python "${EXECUTOR}" \
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
