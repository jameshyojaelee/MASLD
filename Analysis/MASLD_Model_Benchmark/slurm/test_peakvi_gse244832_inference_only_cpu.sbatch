#!/usr/bin/env bash
#SBATCH --job-name=model-check-340
#SBATCH --partition=cpu
#SBATCH --account=nslab
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
OUTPUT=${ROOT}/executions/model-check-340-${SLURM_JOB_ID}
STAGE=${OUTPUT}.staging

if [[ -e ${OUTPUT} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite GSE244832 PeakVI inference-only fixture\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
PYTHONPATH=src:. "${PYTHON}" -m unittest -v \
  tests.unit.test_peakvi_gse244832_inference_only \
  tests.unit.test_commit_peakvi_gse244832_atac_transport_predictions \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
sha256sum \
  scripts/peakvi_gse244832_inference_only.py \
  tests/unit/test_peakvi_gse244832_inference_only.py \
  scripts/commit_peakvi_gse244832_atac_transport_predictions.py \
  tests/unit/test_commit_peakvi_gse244832_atac_transport_predictions.py \
  scripts/commit_gse244832_atac_transport_predictions.py \
  tests/unit/test_commit_gse244832_atac_transport_predictions.py \
  slurm/test_peakvi_gse244832_inference_only_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH=src:. "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]),
    {
        "artifact_class": "peakvi_gse244832_inference_only_fixture",
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse244832",
        "tests": 8,
        "fit_or_refit_performed": False,
        "target_adaptation_performed": False,
        "eligible_donor_lineage_units": 48,
        "missing_predictions_are_nan": True,
        "prediction_commit_firewall_passed": True,
        "rna_assay_read": False,
        "cross_assay_join_consumed": False,
        "outcomes_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${OUTPUT}"
printf 'PeakVI GSE244832 inference-only fixture: %s\n' "${OUTPUT}"
sha256sum "${OUTPUT}/ARTIFACTS.json"
