#!/bin/bash
#SBATCH --job-name=model-check-260
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
CONTRACT="${ROOT}/config/borzoi_grelu_converted_port_weights_only_probe.json"
PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python3.11
TARGET="${ROOT}/executions/model-check-260-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
cd "${ROOT}"
mkdir -p logs executions
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite converted-port weights-only probe: %s\n' \
    "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

PYTHONPATH="src:." "${PYTHON}" -m unittest -v \
  tests.unit.test_probe_borzoi_grelu_converted_port_weights_only \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="src:." "${PYTHON}" \
  scripts/probe_borzoi_grelu_converted_port_weights_only.py \
  --contract "${CONTRACT}" \
  --project-root "${ROOT}" \
  --output "${STAGE}/probe"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONTRACT}" \
  "${ROOT}/scripts/probe_borzoi_grelu_converted_port_weights_only.py" \
  "${ROOT}/tests/unit/test_probe_borzoi_grelu_converted_port_weights_only.py" \
  "${ROOT}/slurm/probe_borzoi_grelu_converted_port_weights_only_cpu.sbatch" \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/weights_only_probe.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "borzoi_grelu_converted_port_weights_only_probe",
        "model_id": receipt["model_id"],
        "checkpoint_members": len(receipt["members"]),
        "replicate_schemas_identical": receipt["replicate_schemas_identical"],
        "strict_local_port_restoration_passed": receipt[
            "strict_local_port_restoration_passed"
        ],
        "unrestricted_torch_load_executed": receipt[
            "unrestricted_torch_load_executed"
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
printf 'Converted-port weights-only probe artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
