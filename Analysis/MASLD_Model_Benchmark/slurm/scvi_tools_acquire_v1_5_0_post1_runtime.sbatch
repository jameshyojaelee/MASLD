#!/usr/bin/env bash
#SBATCH --job-name=masld-scvi150-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
EXPECTED_SCVI_WHEEL_SHA256=d868c86332f819b8d79846ea0e36ec6eef26c3fafeb03ce7a239f1ddefb1c976

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed scvi-tools stage because target exists: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scvi-tools stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scvi-tools runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/downloads" \
  "${STAGE}/validation" \
  "${STAGE}/validation/xdg-cache"

module load python/3.12.3-GCCcore-13.3.0
python -m venv "${STAGE}/venv"
PYTHON="${STAGE}/venv/bin/python"
export PIP_NO_CACHE_DIR=1
export XDG_CACHE_HOME="${STAGE}/validation/xdg-cache"
"${PYTHON}" -m pip install \
  --disable-pip-version-check \
  --upgrade \
  pip==25.2 \
  setuptools==80.9.0 \
  wheel==0.45.1
"${PYTHON}" -m pip download \
  --disable-pip-version-check \
  --no-deps \
  --only-binary=:all: \
  --dest "${STAGE}/downloads" \
  scvi-tools==1.5.0.post1
SCVI_WHEELS=("${STAGE}"/downloads/scvi_tools-1.5.0.post1-*.whl)
if [[ "${#SCVI_WHEELS[@]}" -ne 1 || ! -f "${SCVI_WHEELS[0]}" ]]; then
  printf 'Expected exactly one scvi-tools 1.5.0.post1 wheel\n' >&2
  exit 1
fi
printf '%s  %s\n' "${EXPECTED_SCVI_WHEEL_SHA256}" "${SCVI_WHEELS[0]}" \
  | sha256sum --check --strict
"${PYTHON}" -m pip install \
  --disable-pip-version-check \
  --report "${STAGE}/validation/pip-install-report.json" \
  "${SCVI_WHEELS[0]}"
"${PYTHON}" -m pip check > "${STAGE}/validation/pip-check.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" -m pip list --format=json > "${STAGE}/validation/pip-list.json"
"${PYTHON}" - "${STAGE}" "${EXPECTED_SCVI_WHEEL_SHA256}" <<'PY' \
  > "${STAGE}/validation/import-probe.json"
from pathlib import Path
import json
import platform
import sys

import anndata
import mudata
import numpy
import scvi
import torch

root = Path(sys.argv[1]).resolve(strict=True)
expected = sys.argv[2]
if scvi.__version__ != "1.5.0.post1" or sys.version_info[:2] != (3, 12):
    raise SystemExit("scvi-tools or Python version differs")
from scvi.model import MULTIVI, PEAKVI

payload = {
    "schema_version": "masld-bench-scvi-tools-runtime-probe-v1",
    "status": "pass",
    "python": platform.python_version(),
    "scvi_tools": scvi.__version__,
    "torch": torch.__version__,
    "torch_cuda_build": torch.version.cuda,
    "cuda_visible_on_io_node": torch.cuda.is_available(),
    "numpy": numpy.__version__,
    "anndata": anndata.__version__,
    "mudata": mudata.__version__,
    "multivi_class": MULTIVI.__name__,
    "peakvi_class": PEAKVI.__name__,
    "scvi_wheel_sha256": expected,
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
    "forward_executed": False,
    "gpu_probe_pending": True,
}
print(json.dumps(payload, sort_keys=True))
PY
sha256sum \
  "${ROOT}/slurm/scvi_tools_acquire_v1_5_0_post1_runtime.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${EXPECTED_SCVI_WHEEL_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation/import-probe.json").read_text())
if (
    probe["status"] != "pass"
    or probe["python"].split(".")[:2] != ["3", "12"]
    or probe["scvi_tools"] != "1.5.0.post1"
    or probe["scvi_wheel_sha256"] != sys.argv[2]
    or probe["project_data_read"]
    or probe["outcomes_read"]
    or probe["model_fit_executed"]
    or probe["forward_executed"]
):
    raise SystemExit("scvi-tools runtime admission differs")
freeze_tree(
    root,
    {
        "artifact_class": "scvi_tools_1_5_0_post1_runtime",
        "model_roster": ["peakvi", "multivi"],
        "scvi_tools_version": "1.5.0.post1",
        "scvi_tools_wheel_sha256": sys.argv[2],
        "code_license": "BSD-3-Clause",
        "project_data_read": False,
        "outcomes_read": False,
        "gpu_probe_pending": True,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scvi-tools 1.5.0.post1 runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
