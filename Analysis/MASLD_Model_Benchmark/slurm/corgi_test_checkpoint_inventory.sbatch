#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-schema-test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=executions/masld-corgi-schema-test-%j.out
#SBATCH --error=executions/masld-corgi-schema-test-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/corgi-checkpoint-schema-test-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
SCRIPT_SHA256=2e92448584c148b7e8893245c30e2ffee6c956ca571c908fb4b7279251ba314f
PYTHON=/gpfs/commons/home/jameslee/scImmuneATLAS/.micromamba/envs/immune-atlas/bin/python

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Corgi checkpoint-schema test: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${SCRIPT_SHA256}" \
  "${ROOT}/scripts/corgi_inventory_checkpoint.py" | sha256sum --check --strict

module purge
"${PYTHON}" - "${STAGE}/synthetic.pt" <<'PY'
from pathlib import Path
import sys
import torch

torch.save({
    "state_dict": {
        "aux_encoder.0.weight": torch.zeros((4, 26, 1), dtype=torch.float32),
        "head.weight": torch.ones((2, 4), dtype=torch.float16),
    },
    "epoch": 1,
}, Path(sys.argv[1]))
PY
SIZE="$(stat -c '%s' "${STAGE}/synthetic.pt")"
MD5="$(md5sum "${STAGE}/synthetic.pt" | cut -d' ' -f1)"
SHA256="$(sha256sum "${STAGE}/synthetic.pt" | cut -d' ' -f1)"
"${PYTHON}" "${ROOT}/scripts/corgi_inventory_checkpoint.py" \
  --checkpoint "${STAGE}/synthetic.pt" \
  --expected-size "${SIZE}" \
  --expected-md5 "${MD5}" \
  --expected-sha256 "${SHA256}" \
  --output "${STAGE}/schema.json"
"${PYTHON}" - "${STAGE}/schema.json" <<'PY'
from pathlib import Path
import json
import sys

schema = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
if (
    schema["state_key_count"] != 2
    or schema["wrapper_path"] != ["state_dict"]
    or schema["auxiliary_key_records"] != [{
        "dtype": "torch.float32",
        "key": "aux_encoder.0.weight",
        "numel": 104,
        "shape": [4, 26, 1],
        "storage_bytes": 416,
    }]
    or not schema["safe_loader"]["torch_load_weights_only"]
    or schema["safe_loader"]["corgi_imported"]
):
    raise SystemExit("safe Corgi checkpoint-schema fixture differs")
PY
"${PYTHON}" -m pip freeze > "${STAGE}/pip-freeze.txt"
"${PYTHON}" - "${STAGE}/runtime.json" <<'PY'
from pathlib import Path
import json
import platform
import sys
import torch

Path(sys.argv[1]).write_text(json.dumps({
    "python": platform.python_version(),
    "torch": torch.__version__,
    "torch_cuda": torch.version.cuda,
    "cuda_available": torch.cuda.is_available(),
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]),
    metadata={
        "artifact_class": "corgi_checkpoint_schema_test",
        "status": "pass",
        "checkpoint": "synthetic_only",
        "evaluator_outcomes_exposed": False,
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Corgi checkpoint-schema test: %s\n' "${TARGET}"
