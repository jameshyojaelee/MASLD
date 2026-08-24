#!/usr/bin/env bash
#SBATCH --job-name=masld-scvi150-recover
#SBATCH --partition=cpu
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
FAILED="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-21064550.failed"
FAILED_SOURCE_SHA256=9788a71248410475f18cd60e5169eafbc00ded90d59f0b9407a206e9a109d220
EXPECTED_SCVI_WHEEL_SHA256=d868c86332f819b8d79846ea0e36ec6eef26c3fafeb03ce7a239f1ddefb1c976
TARGET="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed portable scvi-tools stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed portable scvi-tools stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ ! -d "${FAILED}" || -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Portable scvi-tools recovery path contract differs\n' >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/downloads" \
  "${STAGE}/validation"

sha256sum --check --strict "${FAILED}/source.sha256"
SCVI_WHEELS=("${FAILED}"/downloads/scvi_tools-1.5.0.post1-*.whl)
if [[ "${#SCVI_WHEELS[@]}" -ne 1 || ! -f "${SCVI_WHEELS[0]}" ]]; then
  printf 'Expected exactly one recovered scvi-tools wheel\n' >&2
  exit 1
fi
printf '%s  %s\n' "${EXPECTED_SCVI_WHEEL_SHA256}" "${SCVI_WHEELS[0]}" \
  | sha256sum --check --strict
if find "${FAILED}/venv/lib/python3.12/site-packages" -type l -print -quit \
  | grep -q .; then
  printf 'Recovered scvi-tools site-packages contains a symlink\n' >&2
  exit 1
fi
cp -a "${FAILED}/venv/lib/python3.12/site-packages" "${STAGE}/site-packages"
cp -a "${SCVI_WHEELS[0]}" "${STAGE}/downloads/"
cp -a \
  "${FAILED}/validation/import-probe.json" \
  "${FAILED}/validation/pip-check.txt" \
  "${FAILED}/validation/pip-freeze.txt" \
  "${FAILED}/validation/pip-install-report.json" \
  "${FAILED}/validation/pip-list.json" \
  "${STAGE}/validation/"

module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${STAGE}/site-packages" python - "${STAGE}" <<'PY' \
  > "${STAGE}/validation/portable-import-probe.json"
from pathlib import Path
import json
import platform
import sys

import anndata
import mudata
import scvi
import torch
from scvi.model import MULTIVI, PEAKVI

root = Path(sys.argv[1]).resolve(strict=True)
if scvi.__version__ != "1.5.0.post1" or sys.version_info[:2] != (3, 12):
    raise SystemExit("portable scvi-tools or Python version differs")
print(json.dumps({
    "schema_version": "masld-bench-scvi-tools-portable-runtime-probe-v1",
    "status": "pass",
    "python": platform.python_version(),
    "scvi_tools": scvi.__version__,
    "torch": torch.__version__,
    "torch_cuda_build": torch.version.cuda,
    "cuda_visible_on_cpu_node": torch.cuda.is_available(),
    "anndata": anndata.__version__,
    "mudata": mudata.__version__,
    "multivi_class": MULTIVI.__name__,
    "peakvi_class": PEAKVI.__name__,
    "portable_site_packages": True,
    "project_data_read": False,
    "outcomes_read": False,
    "model_fit_executed": False,
    "forward_executed": False,
    "gpu_probe_pending": True,
}, sort_keys=True))
PY
sha256sum \
  "${ROOT}/slurm/scvi_tools_recover_v1_5_0_post1_runtime.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" \
  "${EXPECTED_SCVI_WHEEL_SHA256}" \
  "${FAILED_SOURCE_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation/portable-import-probe.json").read_text())
if (
    probe["status"] != "pass"
    or probe["python"].split(".")[:2] != ["3", "12"]
    or probe["scvi_tools"] != "1.5.0.post1"
    or not probe["portable_site_packages"]
    or probe["project_data_read"]
    or probe["outcomes_read"]
    or probe["model_fit_executed"]
    or probe["forward_executed"]
):
    raise SystemExit("portable scvi-tools runtime admission differs")
freeze_tree(
    root,
    {
        "artifact_class": "scvi_tools_1_5_0_post1_portable_runtime",
        "model_roster": ["peakvi", "multivi"],
        "scvi_tools_version": "1.5.0.post1",
        "scvi_tools_wheel_sha256": sys.argv[2],
        "recovered_failed_campaign_source_sha256": sys.argv[3],
        "failed_campaign_reason": "venv_launcher_symlinks_forbidden_by_immutable_tree",
        "portable_site_packages": True,
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
printf 'Portable scvi-tools 1.5.0.post1 runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
