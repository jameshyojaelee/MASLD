#!/usr/bin/env bash
#SBATCH --job-name=masld-torch260-cu124
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
BASE_SITE="${BASE}/site-packages"
INDEX_URL=https://download.pytorch.org/whl/cu124
TARGET="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed PyTorch CUDA 12.4 stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed PyTorch CUDA 12.4 stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite PyTorch CUDA 12.4 overlay: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/site-packages" \
  "${STAGE}/validation" \
  "${STAGE}/validation/xdg-cache"

printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" \
  | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
export PIP_NO_CACHE_DIR=1
export XDG_CACHE_HOME="${STAGE}/validation/xdg-cache"
PYTHONPATH="${BASE_SITE}" python -m pip install \
  --disable-pip-version-check \
  --index-url "${INDEX_URL}" \
  --report "${STAGE}/validation/pip-install-report.json" \
  --target "${STAGE}/site-packages" \
  torch==2.6.0
if find "${STAGE}/site-packages" -type l -print -quit | grep -q .; then
  printf 'PyTorch CUDA 12.4 overlay contains a symlink\n' >&2
  exit 1
fi
PYTHONPATH="${STAGE}/site-packages:${BASE_SITE}" python - \
  "${BASE_ARTIFACTS_SHA256}" <<'PY' \
  > "${STAGE}/validation/import-probe.json"
import importlib.metadata
import json
import platform
import sys

import scvi
import torch

if (
    sys.version_info[:2] != (3, 12)
    or torch.__version__ != "2.6.0+cu124"
    or torch.version.cuda != "12.4"
    or scvi.__version__ != "1.5.0.post1"
):
    raise SystemExit("PyTorch CUDA 12.4 overlay contract differs")
print(json.dumps({
    "schema_version": "masld-bench-pytorch-cu124-overlay-probe-v1",
    "status": "pass",
    "python": platform.python_version(),
    "torch": torch.__version__,
    "torch_cuda_build": torch.version.cuda,
    "scvi_tools": scvi.__version__,
    "triton": importlib.metadata.version("triton"),
    "cuda_visible_on_io_node": torch.cuda.is_available(),
    "base_runtime_artifacts_sha256": sys.argv[1],
    "official_index_url": "https://download.pytorch.org/whl/cu124",
    "project_data_read": False,
    "outcomes_read": False,
    "forward_executed": False,
    "gpu_probe_pending": True,
}, sort_keys=True))
PY
PYTHONPATH="${STAGE}/site-packages:${BASE_SITE}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/slurm/scvi_tools_prepare_torch2_6_cu124_overlay.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" "${BASE_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/import-probe.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["python"].split(".")[:2] != ["3", "12"]
    or receipt["torch"] != "2.6.0+cu124"
    or receipt["torch_cuda_build"] != "12.4"
    or receipt["scvi_tools"] != "1.5.0.post1"
    or receipt["base_runtime_artifacts_sha256"] != sys.argv[2]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["forward_executed"]
):
    raise SystemExit("PyTorch CUDA 12.4 portable overlay differs")
freeze_tree(
    root,
    {
        "artifact_class": "pytorch_2_6_0_cu124_portable_overlay",
        "base_runtime_artifacts_sha256": sys.argv[2],
        "torch_version": "2.6.0+cu124",
        "torch_cuda_build": "12.4",
        "official_index_url": "https://download.pytorch.org/whl/cu124",
        "project_data_read": False,
        "outcomes_read": False,
        "gpu_probe_pending": True,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'PyTorch 2.6.0 CUDA 12.4 overlay: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
