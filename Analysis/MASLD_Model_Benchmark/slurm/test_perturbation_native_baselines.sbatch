#!/usr/bin/env bash
#SBATCH --job-name=masld-perturb-baseline
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/perturbation-baseline-fixture-${SLURM_JOB_ID}"
if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite baseline fixture: %s\n' "${TARGET}" >&2
  exit 1
fi
mkdir "${TARGET}"
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python -m unittest -v tests.unit.test_perturbation_native_baselines \
  2>&1 | tee "${TARGET}/unittest.log"
sha256sum \
  src/masld_bench/perturbation_baselines.py \
  tests/unit/test_perturbation_native_baselines.py \
  slurm/test_perturbation_native_baselines.sbatch \
  > "${TARGET}/source.sha256"
PYTHONPATH=src python - "${TARGET}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "perturbation_native_baseline_fixture",
    "replicate_first": True,
    "cells_as_replicates": False,
    "systema_source_code_reused": False,
    "gse313774_accessed": False,
}))
PY
sha256sum "${TARGET}/ARTIFACTS.json"
