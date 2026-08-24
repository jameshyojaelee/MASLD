#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-runtime-r1
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
FAILED="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-21064553.staging"
FAILED_IMAGE="${FAILED}/nvidia-pytorch-24.01-py3.sif"
IMAGE_SHA256=eace72d1b831a4834d692134749bb82093c29930fb3a0258ae667b8acb840959
IMAGE_SIZE=10135441408
IMAGE_URI=docker://nvcr.io/nvidia/pytorch:24.01-py3
TARGET="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
IMAGE="${STAGE}/nvidia-pytorch-24.01-py3.sif"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite recovered Corgi runtime artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/runtime" "${STAGE}/validation"
test "$(stat -c '%s' "${FAILED_IMAGE}")" -eq "${IMAGE_SIZE}"
printf '%s  %s\n' "${IMAGE_SHA256}" "${FAILED_IMAGE}" | sha256sum --check --strict
cp -- "${FAILED_IMAGE}" "${IMAGE}"
test "$(stat -c '%s' "${IMAGE}")" -eq "${IMAGE_SIZE}"
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" > "${STAGE}/validation/image.sha256"

module purge
module load singularity/4.2.2
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity inspect --json "${IMAGE}" > "${STAGE}/validation/image.inspect.json"
singularity exec --cleanenv --containall \
  --env PYTHONNOUSERSITE=1 --env HOME=/tmp \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/runtime/pip-freeze.txt"
singularity exec --cleanenv --containall \
  --env PYTHONNOUSERSITE=1 --env HOME=/tmp \
  "${IMAGE}" python - \
  > "${STAGE}/validation/runtime-probe.json" <<'PY'
import importlib.metadata
import json
import platform
import site
import sys

import torch

result = {
    "python": platform.python_version(),
    "python_executable": sys.executable,
    "user_site_enabled": site.ENABLE_USER_SITE,
    "torch": torch.__version__,
    "torch_file": torch.__file__,
    "torch_cuda": torch.version.cuda,
    "cuda_available_during_cpu_acquisition": torch.cuda.is_available(),
    "flash_attn_importable": False,
    "flash_attn_mha_importable": False,
    "status": "pass",
}
if result["user_site_enabled"]:
    raise SystemExit("container user site is unexpectedly enabled")
if "/gpfs/" in str(result["torch_file"]):
    raise SystemExit("host torch leaked into the recovered container")
try:
    import flash_attn

    result["flash_attn_importable"] = True
    result["flash_attn_module_version"] = getattr(flash_attn, "__version__", None)
    try:
        result["flash_attn_distribution_version"] = importlib.metadata.version(
            "flash-attn"
        )
    except importlib.metadata.PackageNotFoundError:
        result["flash_attn_distribution_version"] = None
    try:
        from flash_attn.modules.mha import MHA

        result["flash_attn_mha_importable"] = bool(MHA)
    except Exception as exc:
        result["flash_attn_mha_error"] = f"{type(exc).__name__}: {exc}"
except Exception as exc:
    result["flash_attn_error"] = f"{type(exc).__name__}: {exc}"
print(json.dumps(result, indent=2, sort_keys=True))
PY

cat > "${STAGE}/SOURCE.tsv" <<EOF
field\tvalue
runtime_id\tcorgi_nvidia_pytorch_24_01_py3_recovery
image_uri\t${IMAGE_URI}
image_sha256\t${IMAGE_SHA256}
recovery_source_failed_job\t21064553
recovery_reason\thost_user_site_leaked_without_containall
dependencies_installed_by_campaign\tfalse
checkpoint_loaded\tfalse
outcomes_exposed\tfalse
EOF
sha256sum "${ROOT}/slurm/corgi_recover_ngc_pytorch_24_01_runtime.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${IMAGE_URI}" \
  "${IMAGE_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation/runtime-probe.json").read_text())
if probe.get("status") != "pass":
    raise SystemExit("Corgi recovered runtime inspection did not finish")
if probe.get("user_site_enabled") is not False or "/gpfs/" in str(probe.get("torch_file")):
    raise SystemExit("Corgi recovered runtime is not container-isolated")
image = root / "nvidia-pytorch-24.01-py3.sif"
if sha256_file(image) != sys.argv[3]:
    raise SystemExit("recovered Corgi image differs")
freeze_tree(
    root,
    {
        "artifact_class": "corgi_ngc_pytorch_runtime_acquisition",
        "container_source": sys.argv[2],
        "image_sha256": sys.argv[3],
        "pytorch": probe["torch"],
        "torch_cuda": probe["torch_cuda"],
        "flash_attn_importable": probe["flash_attn_importable"],
        "flash_attn_mha_importable": probe["flash_attn_mha_importable"],
        "container_user_site_enabled": False,
        "dependencies_installed_by_campaign": False,
        "checkpoint_loaded": False,
        "evaluator_outcomes_exposed": False,
        "recovered_from_failed_job": 21064553,
        "status": "pass",
    },
)
PY

mv "${STAGE}" "${TARGET}"
printf 'Recovered Corgi NGC runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
