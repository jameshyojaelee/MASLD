#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-cobolt-test
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:l40s:1
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

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
PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
ENVIRONMENT="${ROOT}/executions/environment-locks/gpu_rna_atac_torch_smoke--c2308fc3b7524081/environment.lock.json"
DATA="${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset/gse296875_rna_atac_smoke_1000.h5"
UPSTREAM="${ROOT}/executions/upstream-acquisitions/cobolt-v1.0.1-21060709/source/epurdom-cobolt-cf5a448/cobolt/model/coboltmodel.py"
OUT="${ROOT}/executions/cobolt-adapter-test-${SLURM_JOB_ID}"

cd "${ROOT}"
printf '%s  %s\n' d4b4a0c60c1388edafbe1c3db91c2370598515083d6d84d1ad23a14df907e4f9 src/masld_bench/adapters/rna_atac_cobolt.py | sha256sum -c -
printf '%s  %s\n' 3627b2d2622d70afddc79ab2b1199e7f90943d97ba5c47dbf5347d8cba05abdd tests/unit/test_rna_atac_cobolt.py | sha256sum -c -
printf '%s  %s\n' 07325273a0df91c23427044175ffd072bce329101f4bf78f042876d31db00d3b tests/scientific/test_rna_atac_cobolt_e2e.py | sha256sum -c -
printf '%s  %s\n' a4f5c3a45ca4069e12ffa05d3efbfcd8a0a636da7515d7d4e9174a914f3f2311 "${UPSTREAM}" | sha256sum -c -

PYTHONPATH=src "${PYTHON}" -m unittest -v tests.unit.test_rna_atac_cobolt
PYTHONPATH=src "${PYTHON}" tests/scientific/test_rna_atac_cobolt_e2e.py \
  --environment-lock "${ENVIRONMENT}" \
  --data "${DATA}" \
  --upstream-source "${UPSTREAM}" \
  --output "${OUT}"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${OUT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "adapter_sha256": "d4b4a0c60c1388edafbe1c3db91c2370598515083d6d84d1ad23a14df907e4f9",
        "artifact_class": "cobolt_adapter_fixture",
        "model_id": "cobolt",
        "result_sha256": sha256_file(root / "result.json"),
        "smoke_only": True,
        "upstream_source_sha256": "a4f5c3a45ca4069e12ffa05d3efbfcd8a0a636da7515d7d4e9174a914f3f2311",
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PY
