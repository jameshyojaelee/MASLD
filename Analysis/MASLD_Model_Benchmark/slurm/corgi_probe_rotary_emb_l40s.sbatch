#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-rotary-gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
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
OVERLAY="${ROOT}/executions/environments/corgi-rotary-emb-0.1-flashattn2.0.4-build-r1-21064845"
OVERLAY_ARTIFACTS_SHA256=c60c4ebf0b32df7d158a3a1d781f5d45fdf4187ab0fba9f4bea35605a92e4c3a
TARGET="${ROOT}/executions/corgi-rotary-emb-l40s-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Corgi rotary L40S probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" \
  "${RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${IMAGE_SHA256}" "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' "${OVERLAY_ARTIFACTS_SHA256}" \
  "${OVERLAY}/ARTIFACTS.json" | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${OVERLAY}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

roots = [Path(value).resolve(strict=True) for value in sys.argv[1:-1]]
payload = {root.name: verify_frozen_tree(root)["metadata"] for root in roots}
Path(sys.argv[-1]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total \
  --format=csv,noheader > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --nv --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:rw" \
  --env PYTHONNOUSERSITE=1 \
  --env "PYTHONPATH=${OVERLAY}/site-packages" \
  "${IMAGE}" python - "${STAGE}/validation/kernel-probe.json" <<'PY'
from pathlib import Path
import json
import site
import sys

import torch
import rotary_emb
from flash_attn.layers.rotary import RotaryEmbedding

if site.ENABLE_USER_SITE:
    raise SystemExit("user site leaked into rotary kernel probe")
if torch.cuda.get_device_name(0) != "NVIDIA L40S":
    raise SystemExit("rotary kernel probe is not on L40S")
embedding = RotaryEmbedding(128, base=20000.0).cuda()
x = torch.randn(
    1, 64, 3, 8, 192,
    device="cuda",
    dtype=torch.bfloat16,
    requires_grad=True,
)
output = embedding(x)
if output.shape != x.shape or not torch.isfinite(output).all():
    raise SystemExit("rotary forward contract differs")
output.float().square().mean().backward()
if x.grad is None or not torch.isfinite(x.grad).all():
    raise SystemExit("rotary backward contract differs")
routes = [
    line.split()
    for line in Path("/proc/net/route").read_text(encoding="utf-8").splitlines()[1:]
    if line.strip()
]
if any(fields and fields[0] != "lo" for fields in routes):
    raise SystemExit("rotary kernel probe has a non-loopback network route")
payload = {
    "schema_version": "masld-bench-corgi-rotary-l40s-v1",
    "status": "pass",
    "torch": torch.__version__,
    "torch_cuda": torch.version.cuda,
    "gpu": torch.cuda.get_device_name(0),
    "gpu_capability": list(torch.cuda.get_device_capability(0)),
    "rotary_emb_file": rotary_emb.__file__,
    "forward": True,
    "backward": True,
    "dtype": "bfloat16",
    "shape": list(output.shape),
    "network_namespace": "none",
    "user_site_enabled": False,
    "evaluator_outcomes_exposed": False,
}
Path(sys.argv[1]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(payload, sort_keys=True))
PY

sha256sum "${ROOT}/slurm/corgi_probe_rotary_emb_l40s.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" \
  "${RUNTIME_ARTIFACTS_SHA256}" "${OVERLAY_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation/kernel-probe.json").read_text())
if (
    probe.get("status") != "pass"
    or probe.get("gpu_capability") != [8, 9]
    or probe.get("forward") is not True
    or probe.get("backward") is not True
):
    raise SystemExit("Corgi rotary L40S kernel probe differs")
freeze_tree(
    root,
    {
        "artifact_class": "corgi_rotary_emb_l40s_probe",
        "runtime_artifacts_sha256": sys.argv[2],
        "overlay_artifacts_sha256": sys.argv[3],
        "gpu": "NVIDIA L40S",
        "gpu_capability": [8, 9],
        "forward_backward_passed": True,
        "network_disabled": True,
        "evaluator_outcomes_exposed": False,
        "status": "pass",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Corgi rotary embedding L40S probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
