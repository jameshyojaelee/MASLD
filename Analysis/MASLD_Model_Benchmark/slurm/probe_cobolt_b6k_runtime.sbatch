#!/usr/bin/env bash
#SBATCH --job-name=model-probe-011
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --gres=gpu:b6k:1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Cobolt B6K runtime feasibility probe.
#
# Re-runs the exact numeric content of the L40S runtime probe (job 21060919) on an
# NVIDIA RTX PRO 6000 Blackwell (sm_120) under torch 2.8.0+cu128, to decide whether
# the Cobolt 25-fit GPU bundle can move off the L40S queue.  Submitted directly
# rather than through the GPU dependency dispatcher: the user authorized B6K in full
# on 2026-08-25, and the dispatcher hardcodes gpu:l40s:1 so it cannot route B6K.
# Job ids from this wrapper never enter the dispatcher roster, so job 21100574 is
# unaffected.
#
# Synthetic tensors only.  Reads no biological data, no development outcome, no
# held-back outcome, and calculates no benchmark metric.  Writes only under executions/,
# which the fixed source record tree walker excludes, so no frozen plan is disturbed.

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
L40S_PROBE="${ROOT}/executions/cobolt-runtime-probe-21060919"
L40S_STATE_SHA256=fa329ed1f4649c1cc375ca042211e09ee74402d4b1ba5df9aaca0a56edbadaeb
ENV_TREE="${ROOT}/executions/environments/scprint2-native-torch2.8.0-v2-21073597"
ENV_ARTIFACTS_SHA256=4474a0d069d09a9a805e7d6e33936af95eeab9fb72d6e3c131e7e7733350842d
PY_BIN="${ENV_TREE}/source/.venv/bin/python"

TARGET="${ROOT}/executions/cobolt-b6k-runtime-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite probe artifact: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}"

cd "${ROOT}"
sha256sum --check --strict <<EOF
${SOURCE_SHA256}  ${SOURCE}
${ACQUISITION_SHA256}  ${ACQUISITION}/ARTIFACTS.json
${L40S_STATE_SHA256}  ${L40S_PROBE}/state_dict.pt
${ENV_ARTIFACTS_SHA256}  ${ENV_TREE}/ARTIFACTS.json
EOF

nvidia-smi > "${STAGE}/nvidia-smi.txt" 2>&1 || true

"${PY_BIN}" - "${SOURCE}" "${STAGE}" "${L40S_PROBE}/state_dict.pt" "${ENV_TREE}" <<'PY'
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import random
import sys
import traceback

import numpy as np
import torch

source = Path(sys.argv[1]).resolve(strict=True)
stage = Path(sys.argv[2]).resolve(strict=True)
l40s_state = Path(sys.argv[3]).resolve(strict=True)
env_tree = sys.argv[4]

result: dict[str, object] = {
    "schema_version": "masld-bench-cobolt-b6k-runtime-probe-v1",
    "source_revision": "cf5a448c6539025346a6393c215ba329cb3a0183",
    "source_tag": "v1.0.1",
}

def fail(stage_name: str, error: BaseException) -> None:
    result["status"] = "fail"
    result["failed_stage"] = stage_name
    result["error"] = f"{type(error).__name__}: {error}"
    (stage / "traceback.txt").write_text(
        "".join(traceback.format_exception(type(error), error, error.__traceback__)),
        encoding="utf-8",
    )
    (stage / "probe.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True, indent=2))
    raise SystemExit(1)

# --- layer guard: prove the intended torch layer is the one that loaded ---
resolved_torch = os.path.realpath(torch.__file__)
result["torch_file"] = resolved_torch
result["torch_version"] = torch.__version__
result["cuda_runtime"] = torch.version.cuda
result["arch_list"] = list(torch.cuda.get_arch_list())
result["layer_guard_passed"] = env_tree.rstrip("/").split("/")[-1] in resolved_torch
if not result["layer_guard_passed"]:
    fail("layer_guard", RuntimeError(f"unexpected torch layer: {resolved_torch}"))
if torch.__version__ != "2.8.0+cu128":
    fail("torch_version", RuntimeError(f"torch is {torch.__version__}"))
if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
    fail("cuda_visibility", RuntimeError("exactly one CUDA device is required"))

device_name = torch.cuda.get_device_name(0)
capability = tuple(int(v) for v in torch.cuda.get_device_capability(0))
result["device_name"] = device_name
result["device_capability"] = list(capability)
result["is_l40s"] = "L40S" in device_name
result["sm_120_native_sass"] = f"sm_{capability[0]}{capability[1]}" in result["arch_list"]
if capability[0] != 12:
    fail("device_identity", RuntimeError(f"expected a Blackwell sm_12x device, got {device_name} {capability}"))
if os.environ.get("CUBLAS_WORKSPACE_CONFIG") != ":4096:8":
    fail("cublas_workspace", RuntimeError("CUBLAS_WORKSPACE_CONFIG=:4096:8 is required"))

# --- determinism switches: the load-bearing unknown on sm_120 ---
try:
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.enabled = False
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
except Exception as error:  # noqa: BLE001
    fail("determinism_switches", error)
result["deterministic_algorithms"] = torch.are_deterministic_algorithms_enabled()
result["cudnn_enabled"] = torch.backends.cudnn.enabled

# --- import the exact upstream v1.0.1 model source ---
try:
    spec = importlib.util.spec_from_file_location("cobolt_v1_0_1_exact_model", source)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import exact Cobolt model source")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
except Exception as error:  # noqa: BLE001
    fail("exact_source_import", error)
result["exact_source_import_passed"] = True


def build_and_train(seed: int):
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
    return model, float(loss.detach().cpu())


# --- check 1: forward/backward equation ---
try:
    model, loss_value = build_and_train(1103)
except Exception as error:  # noqa: BLE001
    fail("forward_backward", error)
result["forward_backward_passed"] = True
result["training_loss"] = loss_value

# in-process determinism: same seed, second construction, identical loss
try:
    _, loss_repeat = build_and_train(1103)
except Exception as error:  # noqa: BLE001
    fail("forward_backward_repeat", error)
result["in_process_seed_repeatability_passed"] = loss_repeat == loss_value
result["training_loss_repeat"] = loss_repeat

# --- check 2: decoder equation parity and query batch invariance ---
try:
    model.eval()
    random.seed(97)
    np.random.seed(97)
    torch.manual_seed(97)
    torch.cuda.manual_seed_all(97)
    query = torch.poisson(torch.full((3, 11), 2.0, device="cuda"))
    with torch.no_grad():
        posterior_mu, _ = model.get_posterior([query, None], [True, False])
        project_profile = torch.softmax(
            torch.softmax(posterior_mu, dim=1) @ model.beta[1], dim=1
        )
        direct_profile = torch.softmax(
            torch.matmul(torch.softmax(posterior_mu, dim=1), model.beta[1]), dim=1
        )
        single_mu, _ = model.get_posterior([query[:1], None], [True, False])
        single_profile = torch.softmax(
            torch.softmax(single_mu, dim=1) @ model.beta[1], dim=1
        )
    if not torch.allclose(project_profile, direct_profile, rtol=0.0, atol=0.0):
        raise RuntimeError("Cobolt decoder equation parity failed")
    result["decoder_equation_parity_passed"] = True
    invariance_error = float(torch.max(torch.abs(project_profile[:1] - single_profile)).cpu())
    result["maximum_query_invariance_error"] = invariance_error
    if not torch.allclose(project_profile[:1], single_profile, rtol=0.0, atol=1e-7):
        raise RuntimeError(
            f"Cobolt frozen query output depends on batch companions (max abs {invariance_error:.3e})"
        )
    result["query_batch_invariance_passed"] = True
    sum_error = float(torch.max(torch.abs(project_profile.sum(dim=1) - 1)).cpu())
    result["maximum_profile_sum_error"] = sum_error
    if sum_error > 1e-6:
        raise RuntimeError("Cobolt decoded profiles do not sum to one")
    result["rna_only_missing_atac_passed"] = True
    result["profile_shape"] = list(project_profile.shape)
except Exception as error:  # noqa: BLE001
    fail("query_invariance", error)

# --- check 3: tensor state restore ---
try:
    state_path = stage / "state_dict.pt"
    torch.save(
        {key: value.detach().cpu() for key, value in model.state_dict().items()},
        state_path,
    )
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
    result["tensor_state_restore_passed"] = True
except Exception as error:  # noqa: BLE001
    fail("tensor_state_restore", error)

# --- check 4: cross-version load of the frozen L40S torch 2.3.1 state dict ---
try:
    legacy = torch.load(l40s_state, map_location="cpu", weights_only=True)
    legacy_model = module.CoboltModel(
        in_channels=[11, 13],
        latent_dim=10,
        n_dataset=[1, 1],
        hidden_dims=[8, 6],
        intercept_adj=False,
        slope_adj=False,
        log=True,
    )
    legacy_model.load_state_dict(legacy, strict=True)
    legacy_model.cuda().eval()
    with torch.no_grad():
        legacy_mu, _ = legacy_model.get_posterior([query, None], [True, False])
        legacy_profile = torch.softmax(
            torch.softmax(legacy_mu, dim=1) @ legacy_model.beta[1], dim=1
        )
    if not torch.isfinite(legacy_profile).all():
        raise RuntimeError("legacy state dict produced non-finite output")
    result["legacy_l40s_state_dict_loadable"] = True
    result["legacy_profile_sum_error"] = float(
        torch.max(torch.abs(legacy_profile.sum(dim=1) - 1)).cpu()
    )
except Exception as error:  # noqa: BLE001
    fail("legacy_state_dict_load", error)

# --- environment shape needed by the project adapter requirements ---
import subprocess

pip = subprocess.run(
    [sys.executable, "-m", "pip", "freeze", "--all"],
    check=False,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
    env={**os.environ, "PYTHONNOUSERSITE": "1"},
)
result["pip_freeze_available"] = pip.returncode == 0
result["pip_freeze_returncode"] = pip.returncode

result["dataset_adjustments_enabled"] = False
result["whole_module_pickle_loaded"] = False
result["status"] = "pass"
(stage / "probe.json").write_text(
    json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(result, sort_keys=True, indent=2))
PY

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${SOURCE_SHA256}" "${ACQUISITION_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "cobolt_b6k_runtime_probe",
        "acquisition_artifacts_sha256": sys.argv[3],
        "model_id": "cobolt",
        "probe_sha256": sha256_file(root / "probe.json"),
        "source_sha256": sys.argv[2],
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Cobolt B6K runtime probe complete: %s\n' "${TARGET}"
