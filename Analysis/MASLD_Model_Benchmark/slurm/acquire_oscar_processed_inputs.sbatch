#!/usr/bin/env bash
#SBATCH --job-name=masld-oscar-inputs
#SBATCH --partition=io
#SBATCH --array=0-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT="${PROJECT}/Analysis/MASLD_Model_Benchmark"
TARGET="${ROOT}/executions/oscar-processed-input-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'OSCAR processed-input acquisition target already exists\n' >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
python scripts/acquire_oscar_processed_inputs.py \
  --task-index "${SLURM_ARRAY_TASK_ID}" \
  --output "${STAGE}"
sha256sum scripts/acquire_oscar_processed_inputs.py slurm/acquire_oscar_processed_inputs.sbatch > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "oscar_author_processed_input_acquisition",
    "source_commit": "57dce224345b3ba25baf585c8a025d4b156d5342",
    "license": "MIT",
    "outcome_model_fit": False,
    "biological_replication_claim": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
