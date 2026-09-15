#!/usr/bin/env bash
#SBATCH --job-name=masld-oscar-topology
#SBATCH --partition=io
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

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT="${PROJECT}/Analysis/MASLD_Model_Benchmark"
SOURCE="${ROOT}/executions/gene-perturbation-source-cra009621_oscar-21066693_3"
TARGET="${ROOT}/executions/oscar-source-topology-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'OSCAR source-topology target already exists\n' >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
python scripts/audit_oscar_source_topology.py \
  --gsa-html "${SOURCE}/raw/CRA009621.html" \
  --output "${STAGE}"
sha256sum scripts/audit_oscar_source_topology.py slurm/audit_oscar_source_topology.sbatch \
  "${SOURCE}/ARTIFACTS.json" > "${STAGE}/source.sha256"
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "oscar_source_topology",
    "article_license": "CC_BY_4.0",
    "production_capture_lanes_as_biological_replicates": False,
    "independent_biological_unit_count": "UNRESOLVED",
    "maximum_independent_units_admitted_for_inference": 1,
    "outcome_model_fit": False,
    "gse313774_accessed": False,
}))
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
