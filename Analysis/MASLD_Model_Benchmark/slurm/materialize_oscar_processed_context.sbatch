#!/usr/bin/env bash
#SBATCH --job-name=masld-oscar-materialize
#SBATCH --partition=io
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT="${PROJECT}/Analysis/MASLD_Model_Benchmark"
PYTHON="${PROJECT}/.envs/perturbation_gears/bin/python"
METADATA="${ROOT}/executions/oscar-processed-input-21068346_10/raw/OSCAR_metadata.csv"
if [[ "${SLURM_ARRAY_TASK_ID}" == "0" ]]; then
  CONTEXT=DM
  CHUNKS=(
    "${ROOT}/executions/oscar-processed-input-21068346_1/raw/OSCAR_DM_expression_matrix.csv.gz.00"
    "${ROOT}/executions/oscar-processed-input-21068346_2/raw/OSCAR_DM_expression_matrix.csv.gz.01"
    "${ROOT}/executions/oscar-processed-input-21068346_3/raw/OSCAR_DM_expression_matrix.csv.gz.02"
    "${ROOT}/executions/oscar-processed-input-21068346_4/raw/OSCAR_DM_expression_matrix.csv.gz.03"
    "${ROOT}/executions/oscar-processed-input-21068346_5/raw/OSCAR_DM_expression_matrix.csv.gz.04"
    "${ROOT}/executions/oscar-processed-input-21068346_6/raw/OSCAR_DM_expression_matrix.csv.gz.05"
  )
else
  CONTEXT=EM
  CHUNKS=(
    "${ROOT}/executions/oscar-processed-input-21068346_7/raw/OSCAR_EM_expression_matrix.csv.gz.00"
    "${ROOT}/executions/oscar-processed-input-21068346_8/raw/OSCAR_EM_expression_matrix.csv.gz.01"
    "${ROOT}/executions/oscar-processed-input-21068346_9/raw/OSCAR_EM_expression_matrix.csv.gz.02"
  )
fi
TARGET="${ROOT}/executions/oscar-${CONTEXT}-materialization-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
STAGE="${TARGET}.staging"
if [[ ! -x "${PYTHON}" || -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'OSCAR materialization precondition differs\n' >&2
  exit 1
fi
cd "${ROOT}"
"${PYTHON}" scripts/materialize_oscar_processed_context.py \
  --context "${CONTEXT}" \
  --metadata "${METADATA}" \
  --chunks "${CHUNKS[@]}" \
  --output "${STAGE}"
sha256sum scripts/materialize_oscar_processed_context.py slurm/materialize_oscar_processed_context.sbatch \
  executions/oscar-processed-input-21068346_*/ARTIFACTS.json > "${STAGE}/source.sha256"
(
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${STAGE}" "${CONTEXT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "oscar_processed_context_materialization",
    "context": sys.argv[2],
    "screen_pool_count": 1,
    "independent_biological_unit_count": "UNRESOLVED",
    "maximum_independent_units_admitted_for_inference": 1,
    "technical_lanes_as_biological_replicates": False,
    "source_outcome_derived_effect_tables_used": False,
    "model_fit_performed": False,
    "gse313774_accessed": False,
}))
PY
)
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
