#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-eval-208
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-eval-208-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-eval-208-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${COBOLT_V8_CANDIDATE:?set COBOLT_V8_CANDIDATE to the exact frozen candidate}
CAMPAIGN_SHA=${COBOLT_V8_CANDIDATE_SHA:?set COBOLT_V8_CANDIDATE_SHA to its ARTIFACTS.json sha256}
SOURCE_LOCK=${COBOLT_V8_EVALUATOR_SOURCE_LOCK:?set COBOLT_V8_EVALUATOR_SOURCE_LOCK to the exact frozen source lock}
SOURCE_LOCK_SHA=${COBOLT_V8_EVALUATOR_SOURCE_LOCK_SHA:?set COBOLT_V8_EVALUATOR_SOURCE_LOCK_SHA to its ARTIFACTS.json sha256}
EVALUATOR_SHA=${COBOLT_V8_EVALUATOR_SHA:?set COBOLT_V8_EVALUATOR_SHA to the v8 entry-point sha256}
EXPECTED_CANDIDATE_PREFIX=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v8--
EXPECTED_SOURCE_PREFIX=${ROOT}/executions/rna-atac-cobolt-5seed-v8-evaluator-source-
for digest in "${CAMPAIGN_SHA}" "${SOURCE_LOCK_SHA}" "${EVALUATOR_SHA}"; do
  if [[ "${digest}" == *[!0-9a-f]* || ${#digest} -ne 64 ]]; then
    printf 'Invalid v8 evaluator digest binding\n' >&2
    exit 1
  fi
done
if [[ "${CANDIDATE}" != "${EXPECTED_CANDIDATE_PREFIX}"* || "${SOURCE_LOCK}" != "${EXPECTED_SOURCE_PREFIX}"* ]]; then
  printf 'Invalid v8 evaluator path binding\n' >&2
  exit 1
fi
EXECUTION_ROOT=${ROOT}/candidates/executions/${CANDIDATE##*/}
EVALUATOR=${ROOT}/src/masld_bench/evaluators/rna_atac_cobolt_5seed_v8_development.py
DATA=${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset/gse296875_rna_atac_smoke_1000.h5
DATA_SHA=a1c05a432cec6894acae20a129ecceb153f191d9911843c86a90e13cd7f8ce14
OUT=${ROOT}/executions/rna-atac-cobolt-5seed-v8-evaluation-${SLURM_JOB_ID}

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${SOURCE_LOCK}" "${CANDIDATE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
from masld_bench.planner import load_frozen_plan

verify_frozen_tree(Path(sys.argv[1]))
candidate = Path(sys.argv[2])
verify_frozen_tree(candidate)
if load_frozen_plan(candidate)["campaign"]["campaign_id"] != "v1-rna-atac-cobolt-5seed-smoke-v8":
    raise SystemExit("candidate campaign is not v8")
PY
printf '%s  %s\n' \
  "${SOURCE_LOCK_SHA}" "${SOURCE_LOCK}/ARTIFACTS.json" \
  "${EVALUATOR_SHA}" "${EVALUATOR}" \
  "${CAMPAIGN_SHA}" "${CANDIDATE}/ARTIFACTS.json" \
  "${DATA_SHA}" "${DATA}" \
  | sha256sum --check --strict

/gpfs/commons/home/jameslee/.local/bin/micromamba run \
  --prefix /gpfs/commons/home/jameslee/micromamba/envs/immune-atlas \
  env \
  PYTHONNOUSERSITE=1 \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  PYTHONPATH=src \
  python -m masld_bench.evaluators.rna_atac_cobolt_5seed_v8_development \
    --candidate "${CANDIDATE}" \
    --execution-root "${EXECUTION_ROOT}" \
    --candidate-sha256 "${CAMPAIGN_SHA}" \
    --h5 "${DATA}" \
    --h5-sha256 "${DATA_SHA}" \
    --output "${OUT}"

install -d -m 0750 "${OUT}/source_snapshot"
install -m 0640 \
  "${SOURCE_LOCK}/source_lock.json" \
  "${OUT}/source_snapshot/source_lock.json"
install -m 0640 "$0" "${OUT}/source_snapshot/evaluation_driver.sbatch"
DRIVER_SHA=$(sha256sum "$0" | cut -d' ' -f1)
PYTHONPATH=src python - \
  "${OUT}" "${CAMPAIGN_SHA}" "${EVALUATOR_SHA}" \
  "${SOURCE_LOCK_SHA}" "${DRIVER_SHA}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

output = Path(sys.argv[1]).resolve(strict=True)
evaluation = output / "evaluation.json"
artifact_hash = freeze_tree(output, {
    "artifact_class": "rna_atac_cobolt_five_seed_development_evaluation",
    "campaign_id": "v1-rna-atac-cobolt-5seed-smoke-v8",
    "campaign_sha256": sys.argv[2],
    "evaluator_sha256": sys.argv[3],
    "evaluator_source_lock_sha256": sys.argv[4],
    "evaluation_driver_sha256": sys.argv[5],
    "evaluation_sha256": sha256_file(evaluation),
    "task_id": "rna_conditioned_atac",
    "smoke_only": True,
    "champion_claim_allowed": False,
})
print(output)
print(artifact_hash)
PY
