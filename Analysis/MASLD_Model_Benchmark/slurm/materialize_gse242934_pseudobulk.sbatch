#!/usr/bin/env bash
#SBATCH --job-name=masld-gse242934-pseudobulk
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
INPUT="${ROOT}/executions/gene-perturbation-training-gse242934_perturbseq-21067795_1"
SOURCE_AUDIT="${ROOT}/executions/gene-perturbation-source-gse242934-21066693_2"
GUIDES="${SOURCE_AUDIT}/raw/GSE242934_Perturbseq_sgrna_lib.txt.gz"
TARGET="${ROOT}/executions/gse242934-pseudobulk-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite GSE242934 materialization: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${INPUT}" "${SOURCE_AUDIT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("source frozen-tree verification passed")
PY
PYTHONPATH=src python scripts/materialize_gse242934_pseudobulk.py \
  --input "${INPUT}" \
  --guide-library "${GUIDES}" \
  --output "${STAGE}"
sha256sum \
  scripts/materialize_gse242934_pseudobulk.py \
  slurm/materialize_gse242934_pseudobulk.sbatch \
  config/datasets/gse242934.toml \
  "${INPUT}/ARTIFACTS.json" \
  "${SOURCE_AUDIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse242934_raw_target_pseudobulk",
    "perturbation_direction": "CRISPRa_activation",
    "biological_unit_count": 1,
    "replicate_inference_allowed": False,
    "normalization_performed": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
