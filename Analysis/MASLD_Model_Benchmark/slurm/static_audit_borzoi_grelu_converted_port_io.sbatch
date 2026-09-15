#!/bin/bash
#SBATCH --job-name=model-data-259
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONTRACT="${ROOT}/config/borzoi_grelu_converted_port_acquisition.json"
FAILED_SOURCE="${ROOT}/executions/model-data-257-21096610.staging/sources"
TARGET="${ROOT}/executions/model-data-259-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
cd "${ROOT}"
mkdir -p logs executions
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite converted-port static audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/sources" "${STAGE}/validation"

for filename in \
  README.md \
  .gitattributes \
  grelu_LICENSE.txt \
  human_state_dict_rep0.h5 \
  human_state_dict_rep1.h5 \
  human_state_dict_rep2.h5 \
  human_state_dict_rep3.h5; do
  cp --reflink=auto --preserve=mode,timestamps \
    "${FAILED_SOURCE}/${filename}" "${STAGE}/sources/${filename}"
done
curl --fail --location --retry 6 --connect-timeout 30 \
  --output "${STAGE}/sources/save_wandb_borzoi_ckpt_human.ipynb.partial" \
  'https://huggingface.co/Genentech/borzoi-model/resolve/9d9e6915f21aba8bfeccca4af2c72deef291af8d/save_wandb_borzoi_ckpt_human.ipynb?download=true'
mv \
  "${STAGE}/sources/save_wandb_borzoi_ckpt_human.ipynb.partial" \
  "${STAGE}/sources/save_wandb_borzoi_ckpt_human.ipynb"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_borzoi_grelu_converted_port \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH="src:." python scripts/audit_borzoi_grelu_converted_port.py \
  --contract "${CONTRACT}" \
  --source-dir "${STAGE}/sources" \
  --output "${STAGE}/audit"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
curl --version > "${STAGE}/validation/curl.version.txt"
sha256sum \
  "${CONTRACT}" \
  "${ROOT}/scripts/audit_borzoi_grelu_converted_port.py" \
  "${ROOT}/tests/unit/test_audit_borzoi_grelu_converted_port.py" \
  "${ROOT}/slurm/static_audit_borzoi_grelu_converted_port_io.sbatch" \
  > "${STAGE}/validation/source.sha256"

PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/converted_port_static_audit.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "borzoi_grelu_converted_port_static_audit",
        "model_id": receipt["model_id"],
        "source_revision": receipt["source_revision"],
        "checkpoint_members": len(receipt["checkpoint_members"]),
        "observed_archive_format": receipt["source_format_incident"][
            "observed_format"
        ],
        "pickle_payload_unpickled": receipt["pickle_payload_unpickled"],
        "unrestricted_torch_load_executed": receipt[
            "unrestricted_torch_load_executed"
        ],
        "torch_weights_only_load_executed": receipt[
            "torch_weights_only_load_executed"
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
printf 'Converted-port static-audit artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
