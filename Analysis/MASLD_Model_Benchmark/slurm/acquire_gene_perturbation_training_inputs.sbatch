#!/usr/bin/env bash
#SBATCH --job-name=masld-perturb-train-inputs
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --array=0-1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%A_%a.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCES=(gse264667_hepg2 gse242934_perturbseq)
SOURCE=${SOURCES[${SLURM_ARRAY_TASK_ID}]}
TARGET="${ROOT}/executions/gene-perturbation-training-${SOURCE}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite training input acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python scripts/acquire_gene_perturbation_training_inputs.py \
  --source "${SOURCE}" \
  --output "${STAGE}"
sha256sum \
  scripts/acquire_gene_perturbation_training_inputs.py \
  slurm/acquire_gene_perturbation_training_inputs.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" "${SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gene_perturbation_training_input_acquisition",
    "source": sys.argv[2],
    "allowlist_exact": True,
    "outcome_model_fit": False,
    "biological_replication_claim": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
