#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)
CANDIDATE_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/spatial_context"

python3 "${SCRIPT_DIR}/01_v1_preservation.py" baseline \
  --project-root "${PROJECT_ROOT}" \
  --candidate-root "${CANDIDATE_ROOT}"
python3 "${SCRIPT_DIR}/02_build_synthetic_contract.py" \
  --project-root "${PROJECT_ROOT}" \
  --candidate-root "${CANDIDATE_ROOT}"
python3 "${SCRIPT_DIR}/03_validate_contract.py" \
  --project-root "${PROJECT_ROOT}" \
  --candidate-root "${CANDIDATE_ROOT}"
python3 "${SCRIPT_DIR}/01_v1_preservation.py" check \
  --project-root "${PROJECT_ROOT}" \
  --candidate-root "${CANDIDATE_ROOT}"

