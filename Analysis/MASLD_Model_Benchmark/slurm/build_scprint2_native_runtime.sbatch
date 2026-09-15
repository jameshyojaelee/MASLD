#!/usr/bin/env bash
#SBATCH --job-name=model-env-302
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
source_artifact="${root}/executions/scprint2-source-21073214"
uv_bin="/gpfs/commons/home/jameslee/.local/bin/uv"
python_bin="/gpfs/commons/home/jameslee/micromamba/envs/scprint2x/bin/python"
target="${root}/executions/environments/scprint2-native-torch2.8.0-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite scPRINT-2 runtime: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

(
  cd "${source_artifact}"
  sha256sum -c SHA256SUMS
)
tar -xf "${source_artifact}/source.tar" -C "${stage}"
printf '%s  %s\n' \
  '088c2bec476ab072b57893d6839dec8014fdddbf8eb1e496559b1ab1338dbd22' \
  "${stage}/source/uv.lock" | sha256sum -c -

export UV_CACHE_DIR="${TMPDIR:-/tmp}/scprint2-uv-cache-${SLURM_JOB_ID}"
export UV_LINK_MODE=copy
export UV_NO_PROGRESS=1
export PYTHONNOUSERSITE=1
"${uv_bin}" sync \
  --frozen \
  --no-dev \
  --python "${python_bin}" \
  --project "${stage}/source"

runtime_python="${stage}/source/.venv/bin/python"
"${runtime_python}" - <<'PY' > "${stage}/runtime_probe.json"
import importlib.metadata
import json
import platform
import sys

import scprint2
import torch

payload = {
    "cuda_available_on_build_node": torch.cuda.is_available(),
    "platform": platform.platform(),
    "python": sys.version,
    "scprint2": importlib.metadata.version("scprint2"),
    "scprint2_module": scprint2.__file__,
    "torch": torch.__version__,
}
print(json.dumps(payload, indent=2, sort_keys=True))
PY
"${uv_bin}" pip freeze --python "${runtime_python}" > "${stage}/pip_freeze.txt"
"${uv_bin}" export \
  --frozen \
  --no-dev \
  --format requirements-txt \
  --project "${stage}/source" \
  --output-file "${stage}/uv_export_requirements.txt"

runtime_probe_sha256="$(sha256sum "${stage}/runtime_probe.json" | cut -d' ' -f1)"
pip_freeze_sha256="$(sha256sum "${stage}/pip_freeze.txt" | cut -d' ' -f1)"
uv_lock_sha256="$(sha256sum "${stage}/source/uv.lock" | cut -d' ' -f1)"
python3 - "${stage}/ARTIFACTS.json" "${runtime_probe_sha256}" \
  "${pip_freeze_sha256}" "${uv_lock_sha256}" <<'PY'
import json
from pathlib import Path
import sys

output, runtime_probe, pip_freeze, uv_lock = sys.argv[1:]
payload = {
    "checkpoint_deserialized": False,
    "runtime_probe_sha256": runtime_probe,
    "pip_freeze_sha256": pip_freeze,
    "schema_version": "masld-bench-scprint2-native-runtime-v1",
    "source_commit": "7ee2de5aaa3492f326f540f703d44ddeaf1f2fd4",
    "source_release": "1.0.3",
    "status": "pass",
    "torch_required": "2.8.0",
    "uv_lock_sha256": uv_lock,
}
Path(output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
(
  cd "${stage}"
  sha256sum ARTIFACTS.json pip_freeze.txt runtime_probe.json uv_export_requirements.txt \
    > SHA256SUMS
)
chmod -R a-w "${stage}"
mv "${stage}" "${target}"
sha256sum "${target}/ARTIFACTS.json"
