#!/usr/bin/env bash
#SBATCH --job-name=masld-gears-core
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT="${PROJECT}/Analysis/MASLD_Model_Benchmark"
PYTHON="${PROJECT}/.envs/perturbation_gears/bin/python"
TARGET="${ROOT}/executions/gears-core-forward-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ ! -x "${PYTHON}" || -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'GEARS core fixture precondition differs\n' >&2
  exit 1
fi
cd "${ROOT}"
"${PYTHON}" scripts/test_gears_core_forward.py --output "${STAGE}"
sha256sum \
  scripts/test_gears_core_forward.py \
  slurm/test_gears_core_forward.sbatch \
  config/artifacts/models/gears/checkpoints.json \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gears_core_forward_fixture",
    "device": "cpu",
    "runtime_network_blocked": True,
    "external_GO_graph_used": False,
    "pickle_or_checkpoint_loaded": False,
    "biological_data_used": False,
    "scientific_admission": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
