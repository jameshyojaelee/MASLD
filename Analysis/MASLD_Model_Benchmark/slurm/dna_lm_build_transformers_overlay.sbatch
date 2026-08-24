#!/usr/bin/env bash
#SBATCH --job-name=masld-dna-lm-runtime
#SBATCH --partition=cpu
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
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
TARGET="${ROOT}/executions/environments/dna-lm-transformers-4.38.1-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite DNA-LM runtime overlay: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/site-packages" \
  "${STAGE}/wheels" \
  "${STAGE}/validation" \
  "${STAGE}/hf-home" \
  "${STAGE}/transformers-cache"
export HF_HOME="${STAGE}/hf-home"
export TRANSFORMERS_CACHE="${STAGE}/transformers-cache"

cat > "${STAGE}/requirements.lock" <<'EOF'
huggingface-hub==0.20.3
regex==2023.12.25
safetensors==0.4.2
tokenizers==0.15.2
transformers==4.38.1
EOF

"${BASE}" -m pip download \
  --disable-pip-version-check \
  --only-binary=:all: \
  --no-deps \
  --dest "${STAGE}/wheels" \
  --requirement "${STAGE}/requirements.lock"
"${BASE}" -m pip install \
  --disable-pip-version-check \
  --no-index \
  --no-deps \
  --target "${STAGE}/site-packages" \
  "${STAGE}"/wheels/*.whl

sha256sum "${STAGE}"/wheels/*.whl > "${STAGE}/validation/wheels.sha256"
PYTHONPATH="${STAGE}/site-packages" "${BASE}" - <<'PY' \
  > "${STAGE}/validation/import_probe.json"
import json
from importlib.metadata import version
import sys

import safetensors
import tokenizers
import transformers

expected = {
    "huggingface-hub": "0.20.3",
    "regex": "2023.12.25",
    "safetensors": "0.4.2",
    "tokenizers": "0.15.2",
    "transformers": "4.38.1",
}
observed = {name: version(name) for name in expected}
if observed != expected:
    raise SystemExit(f"runtime package versions differ: {observed!r}")
print(
    json.dumps(
        {
            "base_python": sys.version,
            "versions": observed,
            "model_code_imported": False,
            "checkpoint_loaded": False,
            "model_forward_executed": False,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
)
PY
PYTHONPATH="${STAGE}/site-packages" "${BASE}" -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
"${BASE}" - <<'PY' > "${STAGE}/validation/base_runtime.json"
import json
import sys
import torch

print(
    json.dumps(
        {
            "python": sys.version,
            "torch": torch.__version__,
            "torch_cuda": torch.version.cuda,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
)
PY
sha256sum \
  "${ROOT}/slurm/dna_lm_build_transformers_overlay.sbatch" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation" / "import_probe.json").read_text())
base = json.loads((root / "validation" / "base_runtime.json").read_text())
if probe["checkpoint_loaded"] or probe["model_forward_executed"]:
    raise SystemExit("runtime admission executed checkpoint or model")
freeze_tree(
    root,
    {
        "artifact_class": "dna_lm_transformers_overlay",
        "base_interpreter": "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python",
        "base_python": base["python"],
        "base_torch": base["torch"],
        "base_torch_cuda": base["torch_cuda"],
        "requirements_lock_sha256": sha256_file(root / "requirements.lock"),
        "runtime_network_allowed": False,
        "runtime_network_used_for_wheel_acquisition_only": True,
    },
)
PY

mv "${STAGE}" "${TARGET}"
printf 'DNA-LM transformers overlay: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
