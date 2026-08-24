#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-rotary-build-r1
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
RUNTIME_ARTIFACTS_SHA256=01de778c17c67d92de9a2e242350e17443dbfe76b7c5a5f8bc447eae5408cd3c
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
IMAGE_SHA256=eace72d1b831a4834d692134749bb82093c29930fb3a0258ae667b8acb840959
FLASH_REV=d30f2e1cd50185c98ed88c0684b4a603f15bee37
PATCH="${ROOT}/patches/flash_attention_v2.0.4_rotary_sm89.patch"
PATCH_SHA256=9646d8a34df8dfc89db163c62ca760389f41964b967cf81967a5f45749afa652
TARGET="${ROOT}/executions/environments/corgi-rotary-emb-0.1-flashattn2.0.4-build-r1-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Corgi rotary overlay: %s\n' "${TARGET}" >&2
  exit 1
fi
if [[ "${PATCH_SHA256}" == TO_BE_FILLED ]]; then
  printf 'Rotary compatibility patch hash is unresolved.\n' >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/downloads" "${STAGE}/source" "${STAGE}/wheel" \
  "${STAGE}/site-packages" "${STAGE}/validation" \
  "${STAGE}/tmp" "${STAGE}/cuda-cache"

printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" \
  "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' "${PATCH_SHA256}" "${PATCH}" | sha256sum --check --strict

curl --fail --location --retry 5 --retry-all-errors --silent --show-error \
  "https://codeload.github.com/Dao-AILab/flash-attention/tar.gz/${FLASH_REV}" \
  --output "${STAGE}/downloads/flash-attention-${FLASH_REV}.tar.gz"
module purge
module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}/downloads/flash-attention-${FLASH_REV}.tar.gz" \
  "${STAGE}/source" "flash-attention-${FLASH_REV}" <<'PY'
from pathlib import Path, PurePosixPath
import sys
import tarfile

archive = Path(sys.argv[1])
destination = Path(sys.argv[2])
expected_root = sys.argv[3]
with tarfile.open(archive, "r:gz") as handle:
    members = handle.getmembers()
    if not members:
        raise SystemExit("empty FlashAttention source archive")
    for member in members:
        path = PurePosixPath(member.name)
        if (
            path.is_absolute()
            or not path.parts
            or path.parts[0] != expected_root
            or ".." in path.parts
            or member.isdev()
            or member.issym()
            or member.islnk()
        ):
            raise SystemExit(f"unsafe FlashAttention source member: {member.name}")
    handle.extractall(destination, members=members)
PY
SOURCE="${STAGE}/source/flash-attention-${FLASH_REV}"
test -f "${SOURCE}/csrc/rotary/setup.py"
test -f "${SOURCE}/csrc/rotary/rotary.cpp"
test -f "${SOURCE}/csrc/rotary/rotary_cuda.cu"
test -f "${SOURCE}/LICENSE"
sha256sum \
  "${SOURCE}/csrc/rotary/setup.py" \
  "${SOURCE}/csrc/rotary/rotary.cpp" \
  "${SOURCE}/csrc/rotary/rotary_cuda.cu" \
  "${SOURCE}/LICENSE" \
  > "${STAGE}/validation/upstream-source.sha256"
patch --batch --forward -d "${SOURCE}" -p1 < "${PATCH}"
rg -n 'compute_89,code=sm_89' "${SOURCE}/csrc/rotary/setup.py" \
  > "${STAGE}/validation/sm89-patch.txt"
sha256sum "${SOURCE}/csrc/rotary/setup.py" \
  > "${STAGE}/validation/patched-source.sha256"

module purge
module load singularity/4.2.2
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --cleanenv --containall \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env MAX_JOBS=2 \
  --env "TMPDIR=${STAGE}/tmp" \
  --env "CUDA_CACHE_PATH=${STAGE}/cuda-cache" \
  "${IMAGE}" bash -lc \
  "cd '${SOURCE}/csrc/rotary' && python -m pip wheel . --no-build-isolation --no-deps --wheel-dir '${STAGE}/wheel'" \
  > "${STAGE}/validation/build.stdout.txt" \
  2> "${STAGE}/validation/build.stderr.txt"
WHEEL="$(find "${STAGE}/wheel" -maxdepth 1 -type f -name 'rotary_emb-0.1-*.whl' -print -quit)"
if [[ -z "${WHEEL}" ]]; then
  printf 'Rotary wheel was not produced.\n' >&2
  exit 1
fi
singularity exec --cleanenv --containall \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip install \
    --no-deps --no-index --target "${STAGE}/site-packages" "${WHEEL}" \
  > "${STAGE}/validation/install.stdout.txt" \
  2> "${STAGE}/validation/install.stderr.txt"

singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env "PYTHONPATH=${STAGE}/site-packages" \
  "${IMAGE}" python - "${STAGE}/validation/import-probe.json" <<'PY'
from pathlib import Path
import json
import site
import sys

import torch
import rotary_emb
from flash_attn.layers.rotary import RotaryEmbedding

if site.ENABLE_USER_SITE:
    raise SystemExit("user site leaked into rotary import probe")
if torch.cuda.is_available():
    raise SystemExit("CPU rotary build unexpectedly exposed a GPU")
if not str(rotary_emb.__file__).endswith(".so"):
    raise SystemExit("rotary overlay did not import a compiled extension")
embedding = RotaryEmbedding(128, base=20000.0)
if embedding is None:
    raise SystemExit("FlashAttention RotaryEmbedding construction failed")
routes = [
    line.split()
    for line in Path("/proc/net/route").read_text(encoding="utf-8").splitlines()[1:]
    if line.strip()
]
if any(fields and fields[0] != "lo" for fields in routes):
    raise SystemExit("rotary kernel probe has a non-loopback network route")
payload = {
    "status": "pass",
    "torch": torch.__version__,
    "torch_cuda": torch.version.cuda,
    "cuda_available": False,
    "rotary_emb_file": rotary_emb.__file__,
    "flash_attn_rotary_import": True,
    "rotary_embedding_constructed": True,
    "gpu_kernel_probe_performed": False,
    "network_namespace": "none",
    "user_site_enabled": False,
}
Path(sys.argv[1]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(payload, sort_keys=True))
PY

sha256sum \
  "${ROOT}/patches/flash_attention_v2.0.4_rotary_sm89.patch" \
  "${ROOT}/slurm/corgi_build_rotary_emb_overlay.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" \
  "${FLASH_REV}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation/import-probe.json").read_text())
if probe.get("status") != "pass" or probe.get("gpu_kernel_probe_performed") is not False:
    raise SystemExit("Corgi rotary overlay import probe differs")
wheel = next((root / "wheel").glob("rotary_emb-0.1-*.whl"))
freeze_tree(
    root,
    {
        "artifact_class": "corgi_rotary_emb_overlay",
        "campaign_revision": 2,
        "supersedes_failed_job": 21064826,
        "runtime_artifacts_sha256": sys.argv[2],
        "flash_attention_revision": sys.argv[3],
        "flash_attention_version": "2.0.4",
        "rotary_emb_version": "0.1",
        "wheel_sha256": sha256_file(wheel),
        "l40s_sm89_build_target_added": True,
        "rotary_kernel_source_modified": False,
        "network_disabled_for_import_probe": True,
        "compiled_extension_import_passed": True,
        "l40s_kernel_probe_performed": False,
        "forward_backward_passed": False,
        "evaluator_outcomes_exposed": False,
        "status": "pass",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Corgi rotary embedding overlay: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
