#!/usr/bin/env bash
#SBATCH --job-name=model-check-337
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
CONFIG="${ROOT}/config/campaigns/peakvi_gse296875_to_gse281367_cell_materialization_revision2_20260825.json"
CONFIG_SHA256=25320873ac05bc0e71d8b3209d30736f7895566f1fbeb3ca0c6b5b068e796900
TARGET="${ROOT}/executions/model-check-337-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite PeakVI cell-hash revision test\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${CONFIG_SHA256}" "${CONFIG}" | sha256sum --check --strict
cd "${ROOT}"
PYTHONPATH="src:." "${PYTHON}" -m unittest -v \
  tests.unit.test_materialize_peakvi_cross_cohort_cells \
  > "${STAGE}/tests.stdout.txt" \
  2> "${STAGE}/tests.stderr.txt"
PYTHONPATH="src:." "${PYTHON}" - "${ROOT}" "${CONFIG}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
from scripts.materialize_peakvi_cross_cohort_cells import validate_config
root = Path(sys.argv[1]).resolve(strict=True)
config = Path(sys.argv[2]).resolve(strict=True)
stage = Path(sys.argv[3]).resolve(strict=True)
validate_config(root, config)
(stage / "receipt.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-peakvi-cell-hash-revision-admission-v1",
            "status": "passed",
            "unit_tests": 4,
            "donor_qualified_cell_hash": True,
            "recurrent_cross_donor_barcode_fixture_passed": True,
            "biological_fragment_values_read": False,
            "condition_values_read": False,
            "outcomes_read": False,
            "metrics_calculated": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
freeze_tree(
    stage,
    {
        "artifact_class": "peakvi_cross_cohort_cell_hash_revision_admission",
        "biological_fragment_values_read": False,
        "outcomes_read": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
