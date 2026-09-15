#!/usr/bin/env bash
#SBATCH --job-name=pool-audit
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
TARGET="${ROOT}/executions/gse281574-pooling-topology-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite pooling audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse281574_pooling_topology
PYTHONPATH="src:." python scripts/audit_gse281574_pooling_topology.py \
  --contract config/evaluation/gse281574_pooling_topology.toml \
  --output "${STAGE}"
sha256sum \
  config/evaluation/gse281574_pooling_topology.toml \
  scripts/audit_gse281574_pooling_topology.py \
  tests/unit/test_audit_gse281574_pooling_topology.py \
  config/evaluation/external_observed_multiome_cohort_expansion.toml \
  config/evaluation/observed_multiome_task.toml \
  OVERALL_PLAN.md \
  slurm/audit_gse281574_pooling_topology_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src:." python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse281574_pooling_topology_audit",
            "outcome_visibility": "public_metadata_only",
            "sealed_outcomes_read": False,
            "biological_matrices_downloaded": False,
            "controlled_data_accessed": False,
            "activates_dataset": False,
            "donor_safe_split_available": False,
            "role": "alcohol_etiology_ood_and_descriptive_only",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE281574 pooling-topology audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
