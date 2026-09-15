#!/usr/bin/env bash
#SBATCH --job-name=model-data-333
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=36:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG=${ROOT}/config/campaigns/peakvi_gse296875_to_gse281367_cell_materialization_20260825.json
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/peakvi-gse296875-source-cells-${SLURM_JOB_ID}
VALID=${ROOT}/executions/peakvi-gse281367-valid-context-cells-${SLURM_JOB_ID}
TEST=${ROOT}/executions/peakvi-gse281367-test-context-cells-${SLURM_JOB_ID}
WORK=${ROOT}/executions/peakvi-cell-materialization-work-${SLURM_JOB_ID}
FAIL=${ROOT}/executions/peakvi-cell-materialization-failed-${SLURM_JOB_ID}

preserve_failure() {
  if [[ -d ${WORK} && ! -e ${FAIL} && ! -e ${WORK}/COMPLETE ]]; then
    mv "${WORK}" "${FAIL}"
  fi
}
trap preserve_failure EXIT

cd "${ROOT}"
PYTHONPATH=src:. "${PYTHON}" -m unittest -v \
  tests.unit.test_materialize_peakvi_cross_cohort_cells
PYTHONPATH=src:. "${PYTHON}" scripts/materialize_peakvi_cross_cohort_cells.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --source-output "${SOURCE}" \
  --target-valid-output "${VALID}" \
  --target-test-output "${TEST}" \
  --work "${WORK}" \
  --workers 16

PYTHONPATH=src:. "${PYTHON}" - "${ROOT}" "${WORK}" "${SOURCE}" "${VALID}" "${TEST}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root, work, *children = map(Path, sys.argv[1:])
for child in children:
    verify_frozen_tree(child)
freeze_tree(
    work,
    {
        "artifact_class": "peakvi_cross_cohort_cell_materialization",
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse281367",
        "native_output_family": "joint_representation",
        "condition_values_read": False,
        "evaluator_outcomes_accessed": False,
        "status": "passed",
    },
)
PY
trap - EXIT
