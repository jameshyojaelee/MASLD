#!/usr/bin/env bash
#SBATCH --job-name=masld-perturb-source
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --array=0-4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCES=(gse238219 gse264667 gse242934 cra009621_oscar igvf_hepatocyte)
SOURCE=${SOURCES[${SLURM_ARRAY_TASK_ID}]}
TARGET="${ROOT}/executions/gene-perturbation-source-${SOURCE}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite source audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python scripts/audit_open_gene_perturbation_source.py \
  --source "${SOURCE}" \
  --output "${STAGE}"
sha256sum \
  scripts/audit_open_gene_perturbation_source.py \
  slurm/audit_open_gene_perturbation_sources.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" "${SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "open_gene_perturbation_source_audit",
    "source": sys.argv[2],
    "outcome_model_fit": False,
    "gse313774_accessed": False,
    "biological_replicate_policy": "source_declared_independent_units_only",
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
