#!/bin/bash
#SBATCH --job-name=model-check-261
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export TOKENIZERS_PARALLELISM=false

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT="${ROOT}/config/borzoi_grelu_to_local_port_mapping_probe.json"
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python3.11
TARGET="${ROOT}/executions/model-check-261-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
REVISION=01b737f2a7da8d251a8912487347c494dcd7bca6
cd "${ROOT}"
mkdir -p logs executions
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite converted-port mapping probe: %s\n' \
    "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/historical_sources" "${STAGE}/validation"

declare -A SOURCE_PATHS=(
  [blocks.py]='src/grelu/model/blocks.py'
  [borzoi.py]='src/grelu/model/trunks/borzoi.py'
  [heads.py]='src/grelu/model/heads.py'
  [layers.py]='src/grelu/model/layers.py'
  [models.py]='src/grelu/model/models.py'
  [position.py]='src/grelu/model/position.py'
)
for filename in $(printf '%s\n' "${!SOURCE_PATHS[@]}" | sort); do
  curl --fail --location --retry 6 --connect-timeout 30 \
    --output "${STAGE}/historical_sources/${filename}.partial" \
    "https://raw.githubusercontent.com/Genentech/gReLU/${REVISION}/${SOURCE_PATHS[$filename]}"
  mv \
    "${STAGE}/historical_sources/${filename}.partial" \
    "${STAGE}/historical_sources/${filename}"
done

PYTHONPATH="src:." "${PYTHON}" -m unittest -v \
  tests.unit.test_reconcile_borzoi_grelu_to_local_port_mapping \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="src:." "${PYTHON}" \
  scripts/reconcile_borzoi_grelu_to_local_port_mapping.py \
  --contract "${CONTRACT}" \
  --project-root "${ROOT}" \
  --historical-sources "${STAGE}/historical_sources" \
  --output "${STAGE}/probe"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONTRACT}" \
  "${ROOT}/scripts/reconcile_borzoi_grelu_to_local_port_mapping.py" \
  "${ROOT}/tests/unit/test_reconcile_borzoi_grelu_to_local_port_mapping.py" \
  "${ROOT}/slurm/reconcile_borzoi_grelu_to_local_port_mapping_cpu.sbatch" \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/mapping_probe.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "borzoi_grelu_to_local_port_semantic_mapping_probe",
        "model_id": receipt["model_id"],
        "checkpoint_members": len(receipt["members"]),
        "mapping_is_bijective": receipt["mapping_is_bijective"],
        "shape_and_dtype_match": receipt["shape_and_dtype_match"],
        "all_four_strict_loads_passed": receipt["all_four_strict_loads_passed"],
        "state_dict_order_used_for_mapping": receipt[
            "state_dict_order_used_for_mapping"
        ],
        "model_forward_executed": receipt["model_forward_executed"],
        "native_borzoi_substitute": receipt["native_borzoi_substitute"],
        "open_champion_eligible": receipt["open_champion_eligible"],
        "terminal_disposition": receipt["terminal_disposition"],
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Converted-port semantic mapping artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
