#!/usr/bin/env bash
#SBATCH --job-name=model-audit-104
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
module load Python/3.11.3-GCCcore-12.3.0

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
export PYTHONPATH="${ROOT}/src"
ADMISSION="${ROOT}/executions/epibert-epcotv2-admission-sources-21064417"
RUNTIME="${ROOT}/executions/environments/epibert-executable-tf2.13-21064522"
SHAPE="${ROOT}/executions/epibert-checkpoint-shape-audit-21065291"
FORWARD="${ROOT}/executions/epibert-released-checkpoint-forward-21065340"
TN5="${ROOT}/executions/gse296875-deduplicated-tn5-bigwigs-21063992"
TARGET="${ROOT}/executions/model-audit-104-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

test ! -e "${TARGET}"
test ! -e "${STAGE}"
mkdir -p "${STAGE}/metadata" "${STAGE}/validation"

curl --fail --location --retry 4 --retry-all-errors \
  --output "${STAGE}/metadata/epcotv2_space_api.json" \
  "https://huggingface.co/api/spaces/luosanj/EPCOTv2/revision/f93a2a944a49797a111f0294709416d38bedbcbe"
curl --fail --location --retry 4 --retry-all-errors \
  --output "${STAGE}/metadata/epcotv2_reference_api.json" \
  "https://huggingface.co/api/datasets/luosanj/epcotv2_data/revision/f445058917800811cb82c38b2c806674cc02b47f"

python3 "${ROOT}/tests/unit/test_audit_epibert_epcotv2_lane.py" \
  > "${STAGE}/validation/unit-test.stdout.txt"
python3 "${ROOT}/scripts/audit_epibert_epcotv2_lane.py" \
  --admission "${ADMISSION}" \
  --runtime "${RUNTIME}" \
  --shape-audit "${SHAPE}" \
  --forward-probe "${FORWARD}" \
  --tn5-bigwigs "${TN5}" \
  --epcot-space-api "${STAGE}/metadata/epcotv2_space_api.json" \
  --epcot-reference-api "${STAGE}/metadata/epcotv2_reference_api.json" \
  --output "${STAGE}/audit"

sha256sum "${ROOT}/scripts/audit_epibert_epcotv2_lane.py" \
  "${ROOT}/tests/unit/test_audit_epibert_epcotv2_lane.py" \
  "${ROOT}/slurm/audit_model_lane_104_cpu.sbatch" \
  > "${STAGE}/source.sha256"
python3 - "${ROOT}" "${STAGE}" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
stage = Path(sys.argv[2])
sys.path.insert(0, str(root / "src"))
from masld_bench.artifacts import freeze_tree

freeze_tree(
    stage,
    metadata={
        "artifact_class": "epibert_epcotv2_lane_audit",
        "models": ["epibert", "epcotv2"],
        "epibert_terminal_disposition": "checkpoint_runtime_ready_native_preprocessing_blocked",
        "epcotv2_terminal_disposition": "blocked_terms_no_execution",
        "gpu_queue_items_created": 0,
        "project_data_read": False,
        "outcomes_read": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
