#!/usr/bin/env bash
#SBATCH --job-name=model-probe-012
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --gres=gpu:b6k:1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Re-validate the three checks named in the Cobolt environment_status string
# -- forward/backward equation, query invariance, tensor state restore -- under
# the NEW cu128 runtime built by slurm/build_cobolt_b6k_runtime.sbatch, on B6K.
#
# The L40S validation (job 21060919, torch 2.3.1+cu121, sm_89) does NOT transfer
# across a torch change or an architecture change and is not reused here.  This
# job re-establishes all three from scratch against the exact upstream v1.0.1
# source, and additionally proves the runtime satisfies every condition
# base._validate_environment enforces at run time -- which is where the scprint2
# cu128 venv fails, since it has no pip.
#
# Submitted directly rather than through the GPU dependency dispatcher: the user
# authorized B6K in full on 2026-08-25, and the dispatcher hardcodes gpu:l40s:1
# so it cannot route B6K.  Job ids from this wrapper never enter the dispatcher
# roster, so job 21100574 is unaffected.
#
# Synthetic tensors only.  Reads no biological data, no development outcome, no
# held-back outcome, and calculates no benchmark metric.
#
# Required exports: COBOLT_B6K_RUNTIME, COBOLT_B6K_RUNTIME_SHA256

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

: "${COBOLT_B6K_RUNTIME:?frozen Cobolt B6K runtime tree is required}"
: "${COBOLT_B6K_RUNTIME_SHA256:?runtime ARTIFACTS.json SHA-256 is required}"

TARGET="${ROOT}/executions/cobolt-b6k-runtime-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite validation artifact: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}"

cd "${ROOT}"
sha256sum --check --strict <<EOF
${SOURCE_SHA256}  ${SOURCE}
${ACQUISITION_SHA256}  ${ACQUISITION}/ARTIFACTS.json
${COBOLT_B6K_RUNTIME_SHA256}  ${COBOLT_B6K_RUNTIME}/ARTIFACTS.json
EOF

PY_BIN="${COBOLT_B6K_RUNTIME}/env/bin/python"
LOCK_DIR=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['environment_lock_directory'])" "${COBOLT_B6K_RUNTIME}/ARTIFACTS.json")
LOCK_JSON="${LOCK_DIR}/environment.lock.json"
printf 'runtime  : %s\nlock dir : %s\n' "${COBOLT_B6K_RUNTIME}" "${LOCK_DIR}"

nvidia-smi > "${STAGE}/nvidia-smi.txt" 2>&1 || true

# Cross-process determinism: the same seeded training run in a separate process.
"${PY_BIN}" - "${SOURCE}" > "${STAGE}/independent_process_loss.txt" <<'PYA'
import importlib.util, random, sys
import numpy as np
import torch

spec = importlib.util.spec_from_file_location("cobolt_v1_0_1_exact_model", sys.argv[1])
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)

torch.use_deterministic_algorithms(True)
torch.backends.cudnn.enabled = False
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True

seed = 1103
random.seed(seed); np.random.seed(seed)
torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
model = module.CoboltModel(
    in_channels=[11, 13], latent_dim=10, n_dataset=[1, 1], hidden_dims=[8, 6],
    intercept_adj=False, slope_adj=False, log=True,
).cuda()
rna = torch.poisson(torch.full((4, 11), 2.0, device="cuda"))
atac = torch.poisson(torch.full((4, 13), 1.0, device="cuda"))
dataset = [torch.zeros(4, dtype=torch.long, device="cuda") for _ in range(2)]
model.train()
latent, recon = model(([rna, atac], dataset), elbo_combn=[[True, False], [False, True], [True, True]])
print(repr(float((latent + recon).detach().cpu())))
PYA

"${PY_BIN}" - "${SOURCE}" "${STAGE}" "${LOCK_JSON}" "${COBOLT_B6K_RUNTIME}" <<'PYB'
from __future__ import annotations

import importlib.metadata
import importlib.util
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import sys
import traceback

import numpy as np
import torch

source = Path(sys.argv[1]).resolve(strict=True)
stage = Path(sys.argv[2]).resolve(strict=True)
lock_json = Path(sys.argv[3]).resolve(strict=True)
runtime_tree = sys.argv[4]

result: dict[str, object] = {
    "schema_version": "masld-bench-cobolt-b6k-runtime-validation-v1",
    "supersedes_status_string": (
        "VALIDATED_JOB_21060919_L40S_TORCH_2_3_1_CUDA_12_1_FORWARD_BACKWARD_"
        "EQUATION_QUERY_INVARIANCE_AND_TENSOR_STATE_RESTORE"
    ),
    "l40s_validation_reused": False,
    "source_revision": "cf5a448c6539025346a6393c215ba329cb3a0183",
    "source_tag": "v1.0.1",
    "runtime_tree": runtime_tree,
}


def fail(stage_name: str, error: BaseException) -> None:
    result["status"] = "fail"
    result["failed_stage"] = stage_name
    result["error"] = f"{type(error).__name__}: {error}"
    (stage / "traceback.txt").write_text(
        "".join(traceback.format_exception(type(error), error, error.__traceback__)),
        encoding="utf-8",
    )
    (stage / "validation.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True, indent=2))
    raise SystemExit(1)


# ---------------------------------------------------------------------------
# Part A: every condition base._validate_environment enforces at run time.
# ---------------------------------------------------------------------------
try:
    lock = json.loads(lock_json.read_text())
    checks: dict[str, bool] = {}
    checks["schema_version"] = lock.get("schema_version") == "masld-bench-python-environment-lock-v1"
    checks["runtime_id"] = lock.get("runtime_id") == "gpu_rna_atac_torch_b6k"
    checks["package_installation_performed_is_false"] = lock.get("package_installation_performed") is False
    checks["pythonnousersite"] = os.environ.get("PYTHONNOUSERSITE") == "1"
    checks["environment_prefix"] = os.path.realpath(sys.prefix) == lock.get("environment_prefix")
    checks["python_executable"] = os.path.realpath(sys.executable) == lock.get("python_executable")
    checks["python_version"] = platform.python_version() == lock.get("python_version")
    critical = lock.get("critical_versions") or {}
    observed_critical = {name: importlib.metadata.version(name) for name in critical}
    checks["critical_versions"] = observed_critical == dict(critical)

    records = sorted(
        {
            f"{str(d.metadata.get('Name')).strip().lower().replace('_', '-')}=={d.version}"
            for d in importlib.metadata.distributions()
            if d.metadata.get("Name")
        }
    )
    distribution_text = "\n".join(records) + "\n"
    checks["python_distributions_sha256"] = (
        sha256(distribution_text.encode()).hexdigest() == lock.get("python_distributions_sha256")
    )

    pip = subprocess.run(
        [sys.executable, "-m", "pip", "freeze", "--all"],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env={**os.environ, "PYTHONNOUSERSITE": "1"},
    )
    checks["pip_freeze_returncode"] = pip.returncode == 0
    checks["pip_freeze_sha256"] = (
        sha256((pip.stdout.rstrip() + "\n").encode()).hexdigest() == lock.get("pip_freeze_sha256")
    )

    result["environment_contract_checks"] = checks
    result["observed_critical_versions"] = observed_critical
    result["environment_contract_passed"] = all(checks.values())
    if not result["environment_contract_passed"]:
        failed = [name for name, ok in checks.items() if not ok]
        raise RuntimeError(f"environment contract failed: {failed}")
except Exception as error:  # noqa: BLE001
    fail("environment_contract", error)

# ---------------------------------------------------------------------------
# Part B: runtime identity.
# ---------------------------------------------------------------------------
try:
    result["torch_version"] = torch.__version__
    result["cuda_runtime"] = torch.version.cuda
    result["arch_list"] = list(torch.cuda.get_arch_list())
    result["torch_file"] = os.path.realpath(torch.__file__)
    result["numpy_version"] = np.__version__
    if runtime_tree.rstrip("/").split("/")[-1] not in result["torch_file"]:
        raise RuntimeError(f"torch did not load from the frozen runtime: {result['torch_file']}")
    if torch.__version__ != "2.8.0+cu128":
        raise RuntimeError(f"torch is {torch.__version__}")
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise RuntimeError("exactly one CUDA device is required")
    device_name = torch.cuda.get_device_name(0)
    capability = tuple(int(v) for v in torch.cuda.get_device_capability(0))
    result["device_name"] = device_name
    result["device_capability"] = list(capability)
    result["sm_120_native_sass"] = f"sm_{capability[0]}{capability[1]}" in result["arch_list"]
    if capability[0] != 12:
        raise RuntimeError(f"expected a Blackwell sm_12x device, got {device_name} {capability}")
    if os.environ.get("CUBLAS_WORKSPACE_CONFIG") != ":4096:8":
        raise RuntimeError("CUBLAS_WORKSPACE_CONFIG=:4096:8 is required")
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.enabled = False
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    result["deterministic_algorithms"] = torch.are_deterministic_algorithms_enabled()
    result["cudnn_enabled"] = torch.backends.cudnn.enabled
except Exception as error:  # noqa: BLE001
    fail("runtime_identity", error)

try:
    spec = importlib.util.spec_from_file_location("cobolt_v1_0_1_exact_model", source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    result["exact_source_import_passed"] = True
except Exception as error:  # noqa: BLE001
    fail("exact_source_import", error)


def build_and_train(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    model = module.CoboltModel(
        in_channels=[11, 13], latent_dim=10, n_dataset=[1, 1], hidden_dims=[8, 6],
        intercept_adj=False, slope_adj=False, log=True,
    ).cuda()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
    rna = torch.poisson(torch.full((4, 11), 2.0, device="cuda"))
    atac = torch.poisson(torch.full((4, 13), 1.0, device="cuda"))
    dataset = [torch.zeros(4, dtype=torch.long, device="cuda") for _ in range(2)]
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


# ---------------------------------------------------------------------------
# Check 1 of 3: forward/backward equation.
# ---------------------------------------------------------------------------
try:
    model, loss_value = build_and_train(1103)
    result["forward_backward_passed"] = True
    result["training_loss"] = loss_value
    independent = float(
        (stage / "independent_process_loss.txt").read_text().strip().strip("'\"")
    )
    result["independent_process_training_loss"] = independent
    result["cross_process_determinism_passed"] = independent == loss_value
    if not result["cross_process_determinism_passed"]:
        raise RuntimeError(
            f"seeded training is not reproducible across processes: {loss_value} vs {independent}"
        )
except Exception as error:  # noqa: BLE001
    fail("forward_backward", error)

# ---------------------------------------------------------------------------
# Check 2 of 3: query invariance (plus decoder equation parity).
# ---------------------------------------------------------------------------
try:
    model.eval()
    random.seed(97); np.random.seed(97)
    torch.manual_seed(97); torch.cuda.manual_seed_all(97)
    query = torch.poisson(torch.full((3, 11), 2.0, device="cuda"))
    with torch.no_grad():
        posterior_mu, _ = model.get_posterior([query, None], [True, False])
        project_profile = torch.softmax(torch.softmax(posterior_mu, dim=1) @ model.beta[1], dim=1)
        direct_profile = torch.softmax(
            torch.matmul(torch.softmax(posterior_mu, dim=1), model.beta[1]), dim=1
        )
        single_mu, _ = model.get_posterior([query[:1], None], [True, False])
        single_profile = torch.softmax(torch.softmax(single_mu, dim=1) @ model.beta[1], dim=1)
    if not torch.allclose(project_profile, direct_profile, rtol=0.0, atol=0.0):
        raise RuntimeError("Cobolt decoder equation parity failed")
    result["decoder_equation_parity_passed"] = True
    invariance_error = float(torch.max(torch.abs(project_profile[:1] - single_profile)).cpu())
    result["maximum_query_invariance_error"] = invariance_error
    if not torch.allclose(project_profile[:1], single_profile, rtol=0.0, atol=1e-7):
        raise RuntimeError(
            f"frozen query output depends on batch companions (max abs {invariance_error:.3e})"
        )
    result["query_batch_invariance_passed"] = True
    sum_error = float(torch.max(torch.abs(project_profile.sum(dim=1) - 1)).cpu())
    result["maximum_profile_sum_error"] = sum_error
    if sum_error > 1e-6:
        raise RuntimeError("decoded profiles do not sum to one")
    result["rna_only_missing_atac_passed"] = True
    result["profile_shape"] = list(project_profile.shape)
except Exception as error:  # noqa: BLE001
    fail("query_invariance", error)

# ---------------------------------------------------------------------------
# Check 3 of 3: tensor state restore.
# ---------------------------------------------------------------------------
try:
    state_path = stage / "state_dict.pt"
    torch.save({k: v.detach().cpu() for k, v in model.state_dict().items()}, state_path)
    restored = module.CoboltModel(
        in_channels=[11, 13], latent_dim=10, n_dataset=[1, 1], hidden_dims=[8, 6],
        intercept_adj=False, slope_adj=False, log=True,
    )
    restored.load_state_dict(torch.load(state_path, map_location="cpu", weights_only=True), strict=True)
    restored.cuda().eval()
    with torch.no_grad():
        restored_mu, _ = restored.get_posterior([query, None], [True, False])
        restored_profile = torch.softmax(
            torch.softmax(restored_mu, dim=1) @ restored.beta[1], dim=1
        )
    if not torch.equal(project_profile, restored_profile):
        raise RuntimeError("tensor-state restore changed predictions")
    result["tensor_state_restore_passed"] = True
except Exception as error:  # noqa: BLE001
    fail("tensor_state_restore", error)

result["dataset_adjustments_enabled"] = False
result["whole_module_pickle_loaded"] = False
result["status"] = "pass"
result["environment_status"] = (
    "VALIDATED_B6K_TORCH_2_8_0_CUDA_12_8_SM120_FORWARD_BACKWARD_EQUATION_"
    "QUERY_INVARIANCE_AND_TENSOR_STATE_RESTORE"
)
(stage / "validation.json").write_text(
    json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(result, sort_keys=True, indent=2))
PYB

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${SOURCE_SHA256}" "${ACQUISITION_SHA256}" "${COBOLT_B6K_RUNTIME_SHA256}" <<'PYC'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "cobolt_b6k_runtime_validation",
        "acquisition_artifacts_sha256": sys.argv[3],
        "model_id": "cobolt",
        "runtime_artifacts_sha256": sys.argv[4],
        "source_sha256": sys.argv[2],
        "validation_sha256": sha256_file(root / "validation.json"),
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PYC

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Cobolt B6K runtime validation complete: %s\n' "${TARGET}"
