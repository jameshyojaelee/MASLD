#!/usr/bin/env bash
#SBATCH --job-name=masld-cobolt-runtime-probe
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --gres=gpu:l40s:1
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ACQUISITION="${ROOT}/executions/upstream-acquisitions/cobolt-v1.0.1-21060709"
SOURCE="${ACQUISITION}/source/epurdom-cobolt-cf5a448/cobolt/model/coboltmodel.py"
SOURCE_SHA256=a4f5c3a45ca4069e12ffa05d3efbfcd8a0a636da7515d7d4e9174a914f3f2311
ACQUISITION_SHA256=73622ff1b36da1da51ab3fc85178f41b5052e1c53e3ae51a75f326bd478b948f
OUT="${ROOT}/executions/cobolt-runtime-probe-${SLURM_JOB_ID}"

mkdir -p "${OUT}"
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}" | sha256sum -c -
printf '%s  %s\n' "${ACQUISITION_SHA256}" "${ACQUISITION}/ARTIFACTS.json" | sha256sum -c -

/usr/bin/env \
  PYTHONNOUSERSITE=1 \
  CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 \
  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python \
  - "${SOURCE}" "${OUT}" <<'PY'
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import random
import sys

import numpy as np
import torch

source = Path(sys.argv[1]).resolve(strict=True)
out = Path(sys.argv[2]).resolve(strict=True)
spec = importlib.util.spec_from_file_location("cobolt_v1_0_1_exact_model", source)
if spec is None or spec.loader is None:
    raise RuntimeError("cannot import exact Cobolt model source")
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)

if torch.__version__ != "2.3.1+cu121":
    raise RuntimeError("torch version differs")
if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
    raise RuntimeError("exactly one CUDA device is required")
device_name = torch.cuda.get_device_name(0)
capability = tuple(torch.cuda.get_device_capability(0))
if "L40S" not in device_name or capability != (8, 9):
    raise RuntimeError("runtime is not the admitted L40S")
torch.use_deterministic_algorithms(True)
torch.backends.cudnn.enabled = False
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True

seed = 1103
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
model = module.CoboltModel(
    in_channels=[11, 13],
    latent_dim=10,
    n_dataset=[1, 1],
    hidden_dims=[8, 6],
    intercept_adj=False,
    slope_adj=False,
    log=True,
).cuda()
optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
rna = torch.poisson(torch.full((4, 11), 2.0, device="cuda"))
atac = torch.poisson(torch.full((4, 13), 1.0, device="cuda"))
dataset = [
    torch.zeros(4, dtype=torch.long, device="cuda"),
    torch.zeros(4, dtype=torch.long, device="cuda"),
]
model.train()
latent_loss, reconstruction_loss = model(
    ([rna, atac], dataset),
    elbo_combn=[[True, False], [False, True], [True, True]],
)
loss = latent_loss + reconstruction_loss
if not torch.isfinite(loss):
    raise RuntimeError("Cobolt training loss is non-finite")
optimizer.zero_grad(set_to_none=True)
loss.backward()
optimizer.step()

model.eval()
query = torch.poisson(torch.full((3, 11), 2.0, device="cuda"))
with torch.no_grad():
    posterior_mu, _ = model.get_posterior([query, None], [True, False])
    project_profile = torch.softmax(
        torch.softmax(posterior_mu, dim=1) @ model.beta[1], dim=1
    )
    direct_profile = torch.softmax(
        torch.matmul(torch.softmax(posterior_mu, dim=1), model.beta[1]),
        dim=1,
    )
    single_mu, _ = model.get_posterior([query[:1], None], [True, False])
    single_profile = torch.softmax(
        torch.softmax(single_mu, dim=1) @ model.beta[1], dim=1
    )
if not torch.allclose(project_profile, direct_profile, rtol=0.0, atol=0.0):
    raise RuntimeError("Cobolt decoder equation parity failed")
if not torch.allclose(project_profile[:1], single_profile, rtol=0.0, atol=1e-7):
    raise RuntimeError("Cobolt frozen query output depends on batch companions")
if not torch.allclose(
    project_profile.sum(dim=1),
    torch.ones(project_profile.shape[0], device="cuda"),
    rtol=0.0,
    atol=1e-6,
):
    raise RuntimeError("Cobolt decoded profiles do not sum to one")

state_path = out / "state_dict.pt"
torch.save({key: value.detach().cpu() for key, value in model.state_dict().items()}, state_path)
restored = module.CoboltModel(
    in_channels=[11, 13],
    latent_dim=10,
    n_dataset=[1, 1],
    hidden_dims=[8, 6],
    intercept_adj=False,
    slope_adj=False,
    log=True,
)
state = torch.load(state_path, map_location="cpu", weights_only=True)
restored.load_state_dict(state, strict=True)
restored.cuda().eval()
with torch.no_grad():
    restored_mu, _ = restored.get_posterior([query, None], [True, False])
    restored_profile = torch.softmax(
        torch.softmax(restored_mu, dim=1) @ restored.beta[1], dim=1
    )
if not torch.equal(project_profile, restored_profile):
    raise RuntimeError("Cobolt tensor-state restore changed predictions")

result = {
    "schema_version": "masld-bench-cobolt-runtime-probe-v1",
    "source_revision": "cf5a448c6539025346a6393c215ba329cb3a0183",
    "source_tag": "v1.0.1",
    "torch_version": torch.__version__,
    "cuda_runtime": torch.version.cuda,
    "device_name": device_name,
    "device_capability": list(capability),
    "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
    "forward_backward_passed": True,
    "rna_only_missing_atac_passed": True,
    "decoder_equation_parity_passed": True,
    "query_batch_invariance_passed": True,
    "tensor_state_restore_passed": True,
    "dataset_adjustments_enabled": False,
    "whole_module_pickle_loaded": False,
    "profile_shape": list(project_profile.shape),
    "maximum_profile_sum_error": float(
        torch.max(torch.abs(project_profile.sum(dim=1) - 1)).cpu()
    ),
}
(out / "probe.json").write_text(
    json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n",
    encoding="utf-8",
)
print(json.dumps(result, sort_keys=True))
PY

PYTHONPATH=src python - "${OUT}" "${SOURCE_SHA256}" "${ACQUISITION_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "cobolt_runtime_probe",
        "acquisition_artifacts_sha256": sys.argv[3],
        "model_id": "cobolt",
        "probe_sha256": sha256_file(root / "probe.json"),
        "source_sha256": sys.argv[2],
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PY
