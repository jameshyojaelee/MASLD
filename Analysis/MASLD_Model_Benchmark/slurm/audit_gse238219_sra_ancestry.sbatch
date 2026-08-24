#!/usr/bin/env bash
#SBATCH --job-name=masld-gse238219-sra
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
EXPERIMENTS=(SRX21157431 SRX21157432 SRX21157433 SRX21157434 SRX21157435)
EXPERIMENT=${EXPERIMENTS[${SLURM_ARRAY_TASK_ID}]}
TARGET="${ROOT}/executions/gse238219-sra-${EXPERIMENT}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite SRA ancestry audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python scripts/audit_gse238219_sra_ancestry.py \
  --experiment "${EXPERIMENT}" \
  --delay-seconds "$((SLURM_ARRAY_TASK_ID * 2))" \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse238219_sra_ancestry.py \
  slurm/audit_gse238219_sra_ancestry.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" "${EXPERIMENT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse238219_sra_ancestry",
    "sra_experiment": sys.argv[2],
    "biological_independence_resolved": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
