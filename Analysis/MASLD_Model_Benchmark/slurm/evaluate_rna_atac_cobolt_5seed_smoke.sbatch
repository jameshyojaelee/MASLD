#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-eval-035
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-eval-035-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-eval-035-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE=${ROOT}/candidates/v1-rna-atac-cobolt-5seed-smoke-v5--162d849f9c6223d5
EXECUTION_ROOT=${ROOT}/candidates/executions/v1-rna-atac-cobolt-5seed-smoke-v5--162d849f9c6223d5
CAMPAIGN_SHA=e7c6afd07c740bc42d3c8ed223be1fdba1b4d72c5eca30adc7c0cf849f4ca17e
SOURCE_LOCK=${ROOT}/executions/rna-atac-cobolt-5seed-evaluator-source-21069723
SOURCE_LOCK_SHA=6165ece8a32c32e2a636eaae488633d43b43b5fadc981d3f7960181d7c587b88
EVALUATOR=${ROOT}/src/masld_bench/evaluators/rna_atac_cobolt_5seed_development.py
EVALUATOR_SHA=43ce6d4d2c08e3191feb6860d73173593586d3febb343b12904f814090fbe8a1
DATA=${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset/gse296875_rna_atac_smoke_1000.h5
DATA_SHA=a1c05a432cec6894acae20a129ecceb153f191d9911843c86a90e13cd7f8ce14
OUT=${ROOT}/executions/rna-atac-cobolt-5seed-evaluation-${SLURM_JOB_ID}

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${SOURCE_LOCK}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
printf '%s  %s\n' \
  "${SOURCE_LOCK_SHA}" "${SOURCE_LOCK}/ARTIFACTS.json" \
  "${EVALUATOR_SHA}" "${EVALUATOR}" \
  "${CAMPAIGN_SHA}" "${CANDIDATE}/ARTIFACTS.json" \
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
  python -m masld_bench.evaluators.rna_atac_cobolt_5seed_development \
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
