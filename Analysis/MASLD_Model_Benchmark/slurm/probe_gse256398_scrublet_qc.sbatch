#!/usr/bin/env bash
#SBATCH --job-name=source-qc-probe-256398
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/gse256398-source-h5-21069407
PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
TARGET=${ROOT}/executions/gse256398-scrublet-qc-probe-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 QC probe: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("GSE256398 source frozen-tree verification passed")
PY
module purge
"${PYTHON}" scripts/probe_gse256398_scrublet_qc.py \
  --h5 "${SOURCE}/h5/GSM8097071_S1_CB_raw_feature_bc_matrix_filtered.h5" \
  --output "${STAGE}/receipt.json" \
  --seed 0
"${PYTHON}" -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum \
  scripts/probe_gse256398_scrublet_qc.py \
  slurm/probe_gse256398_scrublet_qc.sbatch \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt["status"] != "pass"
    or receipt["features_input"] != 36_601
    or receipt["outcomes_or_disease_labels_read"]
    or receipt["model_training_activated"]
):
    raise SystemExit("GSE256398 QC probe receipt differs")
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse256398_scrublet_qc_probe",
            "dataset_id": "gse256398",
            "source_sample_id": "S1",
            "nuclei_input": receipt["nuclei_input"],
            "final_retained": receipt["final_retained"],
            "seed": receipt["seed"],
            "outcome_independent": True,
            "model_training_activated": False,
        },
    )
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE256398 Scrublet/QC probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
