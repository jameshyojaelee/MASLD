#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-eval
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CANDIDATE="${MASLD_CANDIDATE:?frozen candidate path required}"
EXECUTION_ROOT="${MASLD_EXECUTION_ROOT:?scientific execution root required}"
DATA="${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset/gse296875_rna_atac_smoke_1000.h5"
DATA_SHA256=a1c05a432cec6894acae20a129ecceb153f191d9911843c86a90e13cd7f8ce14
OUT="${ROOT}/executions/rna-atac-classical-evaluation-${SLURM_JOB_ID}"

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0

/gpfs/commons/home/jameslee/.local/bin/micromamba run \
  --prefix /gpfs/commons/home/jameslee/micromamba/envs/immune-atlas \
  env \
  PYTHONNOUSERSITE=1 \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  python src/masld_bench/evaluators/rna_atac_development.py \
    --candidate "${CANDIDATE}" \
    --execution-root "${EXECUTION_ROOT}" \
    --h5 "${DATA}" \
    --h5-sha256 "${DATA_SHA256}" \
    --output "${OUT}"

PYTHONPATH=src python - "${OUT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
evaluation = root / "evaluation.json"
freeze_tree(
    root,
    {
        "artifact_class": "rna_atac_development_evaluation",
        "evaluation_sha256": sha256_file(evaluation),
        "task_id": "rna_conditioned_atac",
        "smoke_only": True,
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PY
