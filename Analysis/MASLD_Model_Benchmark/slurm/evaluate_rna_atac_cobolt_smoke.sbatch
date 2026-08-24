#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-cobolt-eval
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
CANDIDATE="${ROOT}/candidates/v1-rna-atac-cobolt-smoke-v2--2f240b2e1fc26f7b"
EXECUTION_ROOT="${ROOT}/candidates/executions/v1-rna-atac-cobolt-smoke-v2--2f240b2e1fc26f7b"
SCPAIR_CANDIDATE="${ROOT}/candidates/v1-rna-atac-scpair-smoke-r2--76f1c71cb509ce1c"
SCPAIR_EXECUTION_ROOT="${ROOT}/candidates/executions/v1-rna-atac-scpair-smoke-r2--76f1c71cb509ce1c"
DATA="${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset/gse296875_rna_atac_smoke_1000.h5"
DATA_SHA256=a1c05a432cec6894acae20a129ecceb153f191d9911843c86a90e13cd7f8ce14
EVALUATOR="${ROOT}/src/masld_bench/evaluators/rna_atac_cobolt_development.py"
EVALUATOR_SHA256=f27486a3eadf57f12cb8644faaa9a8da774fbbf168587613cb9fb3bbc18ca4a1
CAMPAIGN_SHA256=0e2b8880805b9b4b3487b9ab8c7dd2dc384e3be086ed1d25d4293d528297d1bd
SCPAIR_CAMPAIGN_SHA256=9f6b7ae87e9aff96073080ec2bf4d32c423dd1e1400b33f12d4fd11582b195d0
OUT="${ROOT}/executions/rna-atac-cobolt-evaluation-${SLURM_JOB_ID}"

cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
printf '%s  %s\n' "${EVALUATOR_SHA256}" "${EVALUATOR}" | sha256sum -c -
printf '%s  %s\n' "${CAMPAIGN_SHA256}" "${CANDIDATE}/ARTIFACTS.json" | sha256sum -c -
printf '%s  %s\n' "${SCPAIR_CAMPAIGN_SHA256}" "${SCPAIR_CANDIDATE}/ARTIFACTS.json" | sha256sum -c -

mapfile -t attempts < <(find "${EXECUTION_ROOT}/runs" \
  -mindepth 2 -maxdepth 2 -type d -name attempt-001 | sort)
if [[ "${#attempts[@]}" -ne 35 ]]; then
  printf 'Expected 35 Cobolt-campaign attempts; observed %s\n' "${#attempts[@]}" >&2
  exit 1
fi
for attempt in "${attempts[@]}"; do
  PYTHONPATH=src python -m masld_bench.cli run verify --attempt "${attempt}" >/dev/null
done

/gpfs/commons/home/jameslee/.local/bin/micromamba run \
  --prefix /gpfs/commons/home/jameslee/micromamba/envs/immune-atlas \
  env \
  PYTHONNOUSERSITE=1 \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  PYTHONPATH=src \
  python "${EVALUATOR}" \
    --candidate "${CANDIDATE}" \
    --execution-root "${EXECUTION_ROOT}" \
    --scpair-candidate "${SCPAIR_CANDIDATE}" \
    --scpair-execution-root "${SCPAIR_EXECUTION_ROOT}" \
    --h5 "${DATA}" \
    --h5-sha256 "${DATA_SHA256}" \
    --output "${OUT}"

PYTHONPATH=src python - \
  "${OUT}" \
  "${EVALUATOR_SHA256}" \
  "${CAMPAIGN_SHA256}" \
  "${SCPAIR_CAMPAIGN_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
evaluation = root / "evaluation.json"
freeze_tree(
    root,
    {
        "artifact_class": "rna_atac_cobolt_development_evaluation",
        "campaign_sha256": sys.argv[3],
        "scpair_campaign_sha256": sys.argv[4],
        "evaluation_sha256": sha256_file(evaluation),
        "evaluator_sha256": sys.argv[2],
        "task_id": "rna_conditioned_atac",
        "smoke_only": True,
        "champion_claim_allowed": False,
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PY
