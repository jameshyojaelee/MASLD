#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
IMAGE_URI=docker://nvcr.io/nvidia/pytorch:24.01-py3
TARGET="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
IMAGE="${STAGE}/nvidia-pytorch-24.01-py3.sif"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Corgi runtime artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/runtime" "${STAGE}/validation"

module purge
module load singularity/4.2.2
export SINGULARITY_CACHEDIR="${TMPDIR:-/tmp}/singularity-cache-${SLURM_JOB_ID}"
export SINGULARITY_TMPDIR="${TMPDIR:-/tmp}"

singularity pull --disable-cache "${IMAGE}" "${IMAGE_URI}"
test -s "${IMAGE}"
sha256sum "${IMAGE}" > "${STAGE}/validation/image.sha256"
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

import torch

result = {
    "python": platform.python_version(),
    "torch": torch.__version__,
    "torch_cuda": torch.version.cuda,
    "cuda_available_during_cpu_acquisition": torch.cuda.is_available(),
    "flash_attn_importable": False,
    "flash_attn_mha_importable": False,
    "status": "pass",
}
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
runtime_id\tcorgi_nvidia_pytorch_24_01_py3
image_uri\t${IMAGE_URI}
upstream_project\tnvidia/pytorch
upstream_tag\t24.01-py3
intended_models\tcorgi_regular,corgi_plus
dependencies_installed_by_campaign\tfalse
checkpoint_loaded\tfalse
outcomes_exposed\tfalse
EOF
sha256sum "${ROOT}/slurm/corgi_acquire_ngc_pytorch_24_01_runtime.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${IMAGE_URI}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation/runtime-probe.json").read_text())
if probe.get("status") != "pass":
    raise SystemExit("Corgi runtime inspection did not finish")
image = root / "nvidia-pytorch-24.01-py3.sif"
freeze_tree(
    root,
    {
        "artifact_class": "corgi_ngc_pytorch_runtime_acquisition",
        "container_source": sys.argv[2],
        "image_sha256": sha256_file(image),
        "pytorch": probe["torch"],
        "torch_cuda": probe["torch_cuda"],
        "flash_attn_importable": probe["flash_attn_importable"],
        "flash_attn_mha_importable": probe["flash_attn_mha_importable"],
        "dependencies_installed_by_campaign": False,
        "checkpoint_loaded": False,
        "evaluator_outcomes_exposed": False,
        "status": "pass",
    },
)
PY

mv "${STAGE}" "${TARGET}"
printf 'Corgi NGC runtime acquisition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
