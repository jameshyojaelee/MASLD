#!/usr/bin/env bash
#SBATCH --job-name=masld-gse238219-splits
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/gse238219-target-splits-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite target split manifest: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python scripts/build_gse238219_target_split_manifest.py \
  --execution-root "${ROOT}/executions" \
  --materialization-array-job-id 21067095 \
  --output "${STAGE}"
sha256sum \
  scripts/build_gse238219_target_split_manifest.py \
  slurm/build_gse238219_target_split_manifest.sbatch \
  config/tasks/gene_perturbation_response.toml \
  config/splits/gene_perturbation_outer.toml \
  > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gse238219_unseen_target_split_manifest",
    "outcome_independent_role_assignment": True,
    "biological_independence_resolved": False,
    "execution_ready": False,
    "gse281160_nearest_gene_mapping_performed": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
