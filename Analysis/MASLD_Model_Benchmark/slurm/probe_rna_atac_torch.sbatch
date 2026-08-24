#!/usr/bin/env bash
#SBATCH --job-name=masld-rna-atac-torch-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
PYTHON_BIN="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python"
RUN_DIR="${ROOT}/executions/masld-rna-atac-torch-probe-${SLURM_JOB_ID}"

mkdir "${RUN_DIR}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

"${PYTHON_BIN}" - "${RUN_DIR}" <<'PY'
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import sys

import h5py
import numpy
import scipy
import torch

output = Path(sys.argv[1])
torch.use_deterministic_algorithms(True)
torch.backends.cudnn.benchmark = False
torch.manual_seed(20260824)
torch.cuda.manual_seed_all(20260824)

device = torch.device("cuda:0")
model = torch.nn.Sequential(
    torch.nn.Linear(32, 16),
    torch.nn.GELU(),
    torch.nn.Linear(16, 64),
).to(device)
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
features = torch.arange(256, dtype=torch.float32, device=device).reshape(8, 32) / 255.0
targets = torch.arange(512, dtype=torch.float32, device=device).reshape(8, 64) / 511.0

optimizer.zero_grad(set_to_none=True)
prediction = model(features)
loss = torch.nn.functional.mse_loss(prediction, targets)
loss.backward()
optimizer.step()
with torch.no_grad():
    expected = model(features).detach().cpu()

checkpoint = output / "state_dict.pt"
torch.save({key: value.detach().cpu() for key, value in model.state_dict().items()}, checkpoint)
restored = torch.nn.Sequential(
    torch.nn.Linear(32, 16),
    torch.nn.GELU(),
    torch.nn.Linear(16, 64),
).to(device)
state = torch.load(checkpoint, map_location="cpu", weights_only=True)
restored.load_state_dict(state, strict=True)
restored.eval()
with torch.no_grad():
    replay = restored(features).detach().cpu()

if not torch.equal(expected, replay):
    raise RuntimeError("tensor-only checkpoint replay is not bit-identical")
if not torch.isfinite(replay).all():
    raise RuntimeError("probe output is non-finite")

digest = hashlib.sha256(replay.numpy().tobytes()).hexdigest()
record = {
    "schema_version": "masld-bench-rna-atac-torch-probe-v1",
    "status": "pass",
    "python": platform.python_version(),
    "torch": torch.__version__,
    "cuda_runtime": torch.version.cuda,
    "cudnn": torch.backends.cudnn.version(),
    "device_name": torch.cuda.get_device_name(0),
    "device_capability": list(torch.cuda.get_device_capability(0)),
    "numpy": numpy.__version__,
    "scipy": scipy.__version__,
    "h5py": h5py.__version__,
    "loss_before_step": float(loss.detach().cpu()),
    "replay_sha256": digest,
    "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
    "checkpoint_load": "weights_only_true",
}
(output / "probe.json").write_text(
    json.dumps(record, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(record, sort_keys=True))
PY

"${PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/pip_freeze.txt"
sha256sum \
  "${ROOT}/slurm/probe_rna_atac_torch.sbatch" \
  "${RUN_DIR}/probe.json" \
  "${RUN_DIR}/state_dict.pt" \
  > "${RUN_DIR}/sha256.txt"
