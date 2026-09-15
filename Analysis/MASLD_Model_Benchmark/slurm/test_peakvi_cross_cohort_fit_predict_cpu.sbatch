#!/usr/bin/env bash
#SBATCH --job-name=model-check-335
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
RUNTIME="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
TARGET="${ROOT}/executions/model-check-335-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite PeakVI fit-predict test artifact\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
PYTHONPATH="src:." "${RUNTIME}" -m unittest -v \
  tests.unit.test_peakvi_cross_cohort_fit_predict \
  tests.unit.test_commit_gse281367_atac_transport_predictions \
  > "${STAGE}/tests.stdout.txt" \
  2> "${STAGE}/tests.stderr.txt"
sha256sum \
  scripts/peakvi_cross_cohort_fit_predict.py \
  scripts/commit_gse281367_atac_transport_predictions.py \
  tests/unit/test_peakvi_cross_cohort_fit_predict.py \
  tests/unit/test_commit_gse281367_atac_transport_predictions.py \
  slurm/run_peakvi_cross_cohort_bundle.sbatch \
  slurm/test_peakvi_cross_cohort_fit_predict_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src" "${RUNTIME}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
(root / "receipt.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-peakvi-cross-cohort-test-receipt-v1",
            "status": "passed",
            "unit_tests": 8,
            "biological_data_read": False,
            "condition_labels_read": False,
            "development_outcomes_read": False,
            "metrics_calculated": False,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
freeze_tree(
    root,
    {
        "artifact_class": "peakvi_cross_cohort_fit_predict_tests",
        "biological_data_read": False,
        "development_outcomes_read": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
