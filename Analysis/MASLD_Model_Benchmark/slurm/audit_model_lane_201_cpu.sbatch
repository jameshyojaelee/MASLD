#!/usr/bin/env bash
#SBATCH --job-name=model-audit-201
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
module load Python/3.11.3-GCCcore-12.3.0

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
export PYTHONPATH="${ROOT}/src"
CROSS="${ROOT}/executions/scooby-cross-review-21064420"
ARCHIVES="${ROOT}/executions/scooby-released-archive-inventory-21064432"
HEADERS="${ROOT}/executions/scooby-bounded-headers-21064492"
TARGET="${ROOT}/executions/model-audit-201-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

test ! -e "${TARGET}"
test ! -e "${STAGE}"
mkdir -p "${STAGE}/validation"
python3 -m unittest tests.unit.test_audit_scooby_released_lane \
  > "${STAGE}/validation/unit-test.stdout.txt" \
  2> "${STAGE}/validation/unit-test.stderr.txt"
python3 "${ROOT}/scripts/audit_scooby_released_lane.py" \
  --project-root "${ROOT}" \
  --cross-review "${CROSS}" \
  --archive-inventory "${ARCHIVES}" \
  --bounded-headers "${HEADERS}" \
  --output "${STAGE}/audit"
sha256sum \
  "${ROOT}/config/artifacts/models/scooby_neurips/released_variants.json" \
  "${ROOT}/scripts/audit_scooby_released_lane.py" \
  "${ROOT}/tests/unit/test_audit_scooby_released_lane.py" \
  "${ROOT}/slurm/audit_model_lane_201_cpu.sbatch" \
  > "${STAGE}/source.sha256"

python3 - "${ROOT}" "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
stage = Path(sys.argv[2])
sys.path.insert(0, str(root / "src"))
from masld_bench.artifacts import freeze_tree

receipt = json.loads((stage / "audit/released_lane_audit.json").read_text())
freeze_tree(
    stage,
    metadata={
        "artifact_class": "scooby_released_lane_audit",
        "released_checkpoint_count": receipt["released_checkpoint_count"],
        "registered_checkpoint_count": receipt["registered_checkpoint_count"],
        "unregistered_released_variant_count": receipt[
            "unregistered_released_variant_count"
        ],
        "checkpoint_forward_allowed": receipt["checkpoint_forward_allowed"],
        "rna_conditioned_atac_models": receipt["rna_conditioned_atac_models"],
        "sealed_inference_allowed": receipt["sealed_inference_allowed"],
        "checkpoint_payloads_downloaded": False,
        "project_data_read": False,
        "sealed_outcomes_read": False,
        "terminal_disposition": receipt["terminal_disposition"],
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
