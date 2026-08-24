#!/usr/bin/env bash
#SBATCH --job-name=masld-scglue041-runtime
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
SOURCE_ARTIFACT="${ROOT}/executions/upstream-acquisitions/scglue-v0.4.1-21064651"
SOURCE_ARTIFACTS_SHA256=31c64d4da07d433bc7c067595d97340ad90ccaf343423a53c83d932598507282
SOURCE="${SOURCE_ARTIFACT}/source/GLUE-091bdbbdb3dd6d652ed4677200b034ee3c85dae1"
BASE_SITE="${BASE}/site-packages"
TARGET="${ROOT}/executions/environments/scglue-v0.4.1-portable-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed scGLUE runtime stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scGLUE runtime stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scGLUE runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/site-packages" \
  "${STAGE}/wheels" \
  "${STAGE}/validation" \
  "${STAGE}/validation/xdg-cache"

printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE_ARTIFACT}/ARTIFACTS.json" \
  | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${SOURCE_ARTIFACT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
export PIP_NO_CACHE_DIR=1
export XDG_CACHE_HOME="${STAGE}/validation/xdg-cache"
PYTHONPATH="${BASE_SITE}" python -m pip wheel \
  --disable-pip-version-check \
  --no-deps \
  --wheel-dir "${STAGE}/wheels" \
  "${SOURCE}"
SCGLUE_WHEELS=("${STAGE}"/wheels/scglue-0.4.1-*.whl)
if [[ "${#SCGLUE_WHEELS[@]}" -ne 1 || ! -f "${SCGLUE_WHEELS[0]}" ]]; then
  printf 'Expected exactly one locally built scGLUE 0.4.1 wheel\n' >&2
  exit 1
fi
PYTHONPATH="${BASE_SITE}" python -m pip install \
  --disable-pip-version-check \
  --no-deps \
  --target "${STAGE}/site-packages" \
  "${SCGLUE_WHEELS[0]}"
PYTHONPATH="${BASE_SITE}" python -m pip install \
  --disable-pip-version-check \
  --no-deps \
  --report "${STAGE}/validation/dependency-install-report.json" \
  --target "${STAGE}/site-packages" \
  dill \
  igraph \
  leidenalg \
  muon \
  nvidia-ml-py \
  parse \
  pybedtools \
  pysam \
  pytorch-ignite \
  tensorboardX
if find "${STAGE}/site-packages" -type l -print -quit | grep -q .; then
  printf 'scGLUE portable site-packages contains a symlink\n' >&2
  exit 1
fi
PYTHONPATH="${STAGE}/site-packages:${BASE_SITE}" python - \
  "${STAGE}" "${SOURCE_ARTIFACTS_SHA256}" <<'PY' \
  > "${STAGE}/validation/import-probe.json"
from pathlib import Path
from hashlib import sha256
import importlib.metadata
import json
import platform
import sys

import scglue
import torch
from scglue.models.scglue import PairedSCGLUEModel, SCGLUEModel

root = Path(sys.argv[1]).resolve(strict=True)
if scglue.__version__ != "0.4.1" or sys.version_info[:2] != (3, 12):
    raise SystemExit("scGLUE or Python version differs")
packages = [
    "dill", "igraph", "leidenalg", "muon", "nvidia-ml-py", "parse",
    "pybedtools", "pysam", "pytorch-ignite", "tensorboardX",
]
versions = {name: importlib.metadata.version(name) for name in packages}
wheel = next((root / "wheels").glob("scglue-0.4.1-*.whl"))
print(json.dumps({
    "schema_version": "masld-bench-scglue-portable-runtime-probe-v1",
    "status": "pass",
    "python": platform.python_version(),
    "scglue": scglue.__version__,
    "torch": torch.__version__,
    "torch_cuda_build": torch.version.cuda,
    "cuda_visible_on_io_node": torch.cuda.is_available(),
    "scglue_model_class": SCGLUEModel.__name__,
    "paired_scglue_model_class": PairedSCGLUEModel.__name__,
    "resolved_overlay_versions": versions,
    "built_wheel_sha256": sha256(wheel.read_bytes()).hexdigest(),
    "source_artifacts_sha256": sys.argv[2],
    "unsafe_dill_load_called": False,
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
    "forward_executed": False,
    "gpu_probe_pending": True,
}, sort_keys=True))
PY
PYTHONPATH="${STAGE}/site-packages:${BASE_SITE}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/slurm/scglue_prepare_v0_4_1_runtime.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" \
  "${BASE_ARTIFACTS_SHA256}" \
  "${SOURCE_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/import-probe.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["python"].split(".")[:2] != ["3", "12"]
    or receipt["scglue"] != "0.4.1"
    or receipt["source_artifacts_sha256"] != sys.argv[3]
    or receipt["unsafe_dill_load_called"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["model_fit_executed"]
    or receipt["forward_executed"]
):
    raise SystemExit("scGLUE portable runtime contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "scglue_v0_4_1_portable_runtime",
        "base_runtime_artifacts_sha256": sys.argv[2],
        "source_artifacts_sha256": sys.argv[3],
        "scglue_version": "0.4.1",
        "code_license": "MIT",
        "unsafe_dill_load_called": False,
        "project_data_read": False,
        "outcomes_read": False,
        "gpu_probe_pending": True,
        "profile_task_status": "latent_only",
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'Portable scGLUE v0.4.1 runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
