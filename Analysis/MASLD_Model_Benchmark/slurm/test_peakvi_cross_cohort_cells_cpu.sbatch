#!/usr/bin/env bash
#SBATCH --job-name=model-check-334
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
OUTPUT=${ROOT}/executions/model-check-334-${SLURM_JOB_ID}

if [[ -e ${OUTPUT} ]]; then
  printf 'Refusing to overwrite %s\n' "${OUTPUT}" >&2
  exit 1
fi
install -d -m 0750 "${OUTPUT}"
cd "${ROOT}"
PYTHONPATH=src:. "${PYTHON}" -m unittest -v \
  tests.unit.test_materialize_peakvi_cross_cohort_cells \
  >"${OUTPUT}/tests.stdout.txt" 2>"${OUTPUT}/tests.stderr.txt"
PYTHONPATH=src:. "${PYTHON}" - "${ROOT}" "${OUTPUT}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from scripts.materialize_peakvi_cross_cohort_cells import validate_config

root, output = map(Path, sys.argv[1:])
config = root / "config/campaigns/peakvi_gse296875_to_gse281367_cell_materialization_20260825.json"
value, paths = validate_config(root, config)
receipt = {
    "schema_version": "masld-bench-peakvi-cross-cohort-cell-materialization-fixture-v1",
    "status": "passed",
    "campaign_id": value["campaign_id"],
    "authorities_verified": sorted(paths),
    "unit_tests": 3,
    "biological_fragment_values_read": False,
    "condition_values_read": False,
    "evaluator_outcomes_read": False,
    "metrics_calculated": False,
}
(output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
freeze_tree(
    output,
    {
        "artifact_class": "peakvi_cross_cohort_cell_materialization_fixture",
        "biological_fragment_values_read": False,
        "condition_values_read": False,
        "evaluator_outcomes_accessed": False,
        "status": "passed",
    },
)
PY
