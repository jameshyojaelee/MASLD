#!/usr/bin/env bash
#SBATCH --job-name=masld-gse267145-rna-unit
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
MATRIX="${ROOT}/executions/gse267145-matrix-audit-21064333"
JOIN="${ROOT}/executions/gse267145-authoritative-join-21064930"
TARGET="${ROOT}/executions/gse267145-rna-measurement-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite GSE267145 RNA audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_rna_measurement
PYTHONPATH="src" python - "${MATRIX}" "${JOIN}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY
PYTHONPATH="src:." python scripts/audit_gse267145_rna_measurement.py \
  --matrix-audit "${MATRIX}" \
  --matrix-artifacts-sha256 4e0d17f392530f4e6f0b54346d5d9333e450b7606d9698fbf7253a4f2a0e0417 \
  --join-audit "${JOIN}" \
  --join-artifacts-sha256 e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580 \
  --pisces-url https://github.com/Novartis/pisces/archive/refs/tags/v0.1.5.4.tar.gz \
  --pisces-ref-url https://api.github.com/repos/Novartis/pisces/git/ref/tags/v0.1.5.4 \
  --pisces-tag v0.1.5.4 \
  --pisces-commit 9936079ac75d4b75be95bad5bc962465e8c5f458 \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse267145_rna_measurement.py \
  tests/unit/test_audit_gse267145_rna_measurement.py \
  slurm/audit_gse267145_rna_measurement.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse267145_rna_measurement_audit",
            "cohort_family_ids": ["gse267145_znf469_human_liver"],
            "continuous_expression_input_allowed": True,
            "raw_integer_count_likelihood_allowed": False,
            "model_training_activated": False,
            "public_pisces_generator_exact": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 RNA measurement audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
