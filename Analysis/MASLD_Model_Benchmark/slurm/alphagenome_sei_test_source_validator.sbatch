#!/usr/bin/env bash
#SBATCH --job-name=masld-alpha-sei-test
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

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/alphagenome-sei-admission-sources-21064505.staging/sources"
TARGET="${ROOT}/executions/alphagenome-sei-source-validator-test-${SLURM_JOB_ID}"

if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite source-validator test: %s\n' "${TARGET}" >&2
  exit 1
fi
module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/alphagenome_sei_validate_sources.py" \
  --config "${ROOT}/config/alphagenome_sei_admission.json" \
  --source-dir "${SOURCE}" \
  --project-root "${ROOT}" \
  --output "${TARGET}"
python - "${TARGET}/source_admission_receipt.json" <<'PY'
import json
from pathlib import Path
import sys

receipt = json.loads(Path(sys.argv[1]).read_text())
if (
    receipt["status"] != "pass"
    or receipt["checkpoint_bytes_downloaded"]
    or receipt["checkpoint_bytes_deserialized"]
    or receipt["model_forward_executed"]
    or receipt["model_dispositions"]["alphagenome_all_folds"]["next_model_job_allowed"]
    or not receipt["model_dispositions"]["sei"]["safe_archive_acquisition_allowed"]
):
    raise SystemExit("source admission receipt differs")
print("AlphaGenome/Sei source validator passed frozen-source test")
PY
