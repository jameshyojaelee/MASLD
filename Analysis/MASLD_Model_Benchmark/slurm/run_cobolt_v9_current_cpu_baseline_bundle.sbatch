#!/usr/bin/env bash
#SBATCH --job-name=model-training-216
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=18:00:00
#SBATCH --signal=B:USR1@3600
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-216-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-training-216-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
: "${COBOLT_V9_CPU_VALIDATION:?exact CPU validation artifact is required}"
: "${COBOLT_V9_CPU_VALIDATION_SHA:?exact CPU validation ARTIFACTS SHA-256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58
CANDIDATE_SHA=ac24e1a8da4fde1c35aafa1370659155eded36811099b40f3ee39a4fc2866862
STABILITY=${ROOT}/executions/cobolt-v9-frozen-candidate-stability-21087850
STABILITY_SHA=b8b6003576becb93b530c20465d9022ee334bae5bbdcb660ba340836d0831233
EXECUTOR=${ROOT}/scripts/run_cobolt_v9_cpu_baseline_bundle.py
EXECUTOR_SHA=9cca3180174a458a06c713c9b473ad944f26719bed24f10e80ad3675c00d2dc5
EXECUTION_ROOT=${ROOT}/executions/campaign-runs/v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58--cpu-baselines-v2
RECEIPT_ROOT=${ROOT}/executions/cobolt-v9-current-cpu-baseline-bundles

for digest in "${CANDIDATE_SHA}" "${STABILITY_SHA}" "${EXECUTOR_SHA}" "${COBOLT_V9_CPU_VALIDATION_SHA}"; do
  if [[ ${digest} == *[!0-9a-f]* || ${#digest} -ne 64 ]]; then
    printf 'Invalid SHA-256 binding.\n' >&2
    exit 1
  fi
done
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPYCACHEPREFIX=${TMPDIR:-/tmp}/model-training-216-${SLURM_JOB_ID}
cd "${ROOT}"
printf '%s  %s\n' "${CANDIDATE_SHA}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${STABILITY_SHA}" "${STABILITY}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${EXECUTOR_SHA}" "${EXECUTOR}" | sha256sum --check --strict
printf '%s  %s\n' "${COBOLT_V9_CPU_VALIDATION_SHA}" "${COBOLT_V9_CPU_VALIDATION}/ARTIFACTS.json" | sha256sum --check --strict
export PYTHONPATH=${ROOT}:${ROOT}/src
python - "${COBOLT_V9_CPU_VALIDATION}" "${CANDIDATE_SHA}" "${STABILITY_SHA}" "${EXECUTOR_SHA}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
verify_frozen_tree(root)
receipt = json.loads((root / "validation_receipt.json").read_text())
if (
    receipt.get("status") != "pass"
    or receipt.get("candidate_manifest_sha256") != sys.argv[2]
    or receipt.get("stability_manifest_sha256") != sys.argv[3]
    or receipt.get("executor_sha256") != sys.argv[4]
    or receipt.get("preflight_runs") != 150
    or receipt.get("execution_started") is not False
    or receipt.get("outcomes_opened") is not False
):
    raise SystemExit("Cobolt CPU validation binding differs")
PY
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
