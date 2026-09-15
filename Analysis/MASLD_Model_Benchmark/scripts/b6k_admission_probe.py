#!/usr/bin/env python3
"""Outcome-blind GPU inclusion probes for the B6K (sm_120) hardware question.

The campaign runtime is pinned to ``torch 2.6.0+cu124``, whose embedded
``CUDA_ARCH_FLAGS`` are ``sm_50 sm_60 sm_70 sm_75 sm_80 sm_86 sm_90`` with no
``compute_*`` PTX entry.  Whether that build can execute on an sm_120 device is
an empirical question, and probe 1 is the one that answers it.

Every probe uses synthetic tensors.  No biological data, no development
outcome, no held-back outcome, and no benchmark metric is read or computed here.

The same file runs unchanged in all three cells of the inclusion design:

    cu124 / L40S    reference    the campaign's frozen runtime
    cu130 / L40S    control      isolates the torch-version change
    cu130 / B6K     test         isolates the architecture change

The runtime is selected by the caller's ``PYTHONPATH``, never by this script,
so the three cells differ in exactly one variable at a time.  ``--cell-id``
only labels the emitted output file.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import platform
import sys
import traceback
from typing import Any, Callable

import torch
from torch import nn


PROBE_SCHEMA = "masld-bench-b6k-admission-probe-v1"
COMPARE_SCHEMA = "masld-bench-b6k-numeric-tolerance-comparison-v1"

# Probe 4/5 train this many optimizer steps.  Probe 5 splits the run in half
# and must reproduce the uninterrupted trajectory exactly.
TRAIN_STEPS = 24


class ProbeError(RuntimeError):
    """Raised when a probe cannot be executed or its result cannot be trusted."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def tensor_digest(value: torch.Tensor) -> str:
    """Content hash of a tensor, taken in float64 so it is dtype-comparable."""

    flat = value.detach().to("cpu", torch.float64).contiguous().flatten()
    return sha256(flat.numpy().tobytes()).hexdigest()


def _try(label: str, action: Callable[[], Any], record: dict[str, Any]) -> None:
    """Apply an optional backend knob and record whether it actually took.

    torch 2.6 and torch 2.13 do not expose the same precision switches.  A
    silently missing knob would change the numbers the tolerance probe reports,
    so every attempt is recorded rather than swallowed.
    """

    try:
        action()
    except Exception as error:  # noqa: BLE001 - the failure itself is the datum
        record[label] = f"unavailable: {type(error).__name__}"
    else:
        record[label] = "applied"


def pin_determinism(seed: int) -> dict[str, Any]:
    """Pin every knob that would otherwise let the cells differ for free.

    TF32 is the important one.  It is an Ampere-and-later fast path whose
    availability and default differ by architecture and by torch version, and
    leaving it on would make the tolerance probe measure the TF32 policy rather
    than the architecture.
    """

    applied: dict[str, Any] = {}
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    _try(
        "use_deterministic_algorithms",
        lambda: torch.use_deterministic_algorithms(True),
        applied,
    )
    _try(
        "cudnn_deterministic",
        lambda: setattr(torch.backends.cudnn, "deterministic", True),
        applied,
    )
    _try(
        "cudnn_benchmark_off",
        lambda: setattr(torch.backends.cudnn, "benchmark", False),
        applied,
    )
    _try(
        "matmul_tf32_off",
        lambda: setattr(torch.backends.cuda.matmul, "allow_tf32", False),
        applied,
    )
    _try(
        "cudnn_tf32_off",
        lambda: setattr(torch.backends.cudnn, "allow_tf32", False),
        applied,
    )
    _try(
        "fp32_matmul_precision_highest",
        lambda: torch.set_float32_matmul_precision("highest"),
        applied,
    )
    return applied


def arch_flags() -> list[str]:
    """Compile-time arch list, readable without a visible device."""

    try:
        return torch._C._cuda_getArchFlags().split()
    except Exception:  # noqa: BLE001 - a CPU-only build has no flags at all
        return []


def driver_version_from_proc() -> str | None:
    path = Path("/proc/driver/nvidia/version")
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError:
        return None


# --------------------------------------------------------------------------
# Probe 1: import and architecture
# --------------------------------------------------------------------------
def probe_import_architecture() -> dict[str, Any]:
    record: dict[str, Any] = {
        "torch_version": torch.__version__,
        "torch_file": str(Path(torch.__file__).resolve()),
        "torch_cuda_build": torch.version.cuda,
        "torch_arch_list": arch_flags(),
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "cuda_available": bool(torch.cuda.is_available()),
        "nvidia_driver_banner": driver_version_from_proc(),
    }
    try:
        record["cuda_driver_api_version"] = torch._C._cuda_getDriverVersion()
    except Exception as error:  # noqa: BLE001
        record["cuda_driver_api_version"] = f"unavailable: {type(error).__name__}"
    if not record["cuda_available"]:
        record["status"] = "no_cuda_device_visible"
        return record
    capability = torch.cuda.get_device_capability()
    record["device_capability"] = list(capability)
    record["device_sm"] = f"sm_{capability[0]}{capability[1]}"
    record["device_name"] = torch.cuda.get_device_name()
    properties = torch.cuda.get_device_properties(0)
    record["device_total_memory_bytes"] = int(properties.total_memory)
    record["device_multi_processor_count"] = int(properties.multi_processor_count)
    try:
        record["cudnn_version"] = torch.backends.cudnn.version()
    except Exception as error:  # noqa: BLE001
        record["cudnn_version"] = f"unavailable: {type(error).__name__}"
    # The load-bearing question: does this build carry SASS for this device, or
    # PTX it could JIT from?  Neither is fatal on its own; probe 2 decides.
    record["device_sm_in_arch_list"] = record["device_sm"] in record["torch_arch_list"]
    record["arch_list_has_ptx"] = any(
        flag.startswith("compute_") for flag in record["torch_arch_list"]
    )
    record.update(execution_path(capability, record["torch_arch_list"]))
    record["status"] = "pass"
    return record


def parse_arch(flag: str) -> tuple[int, int] | None:
    """Split ``sm_86`` into (8, 6) and ``sm_120`` into (12, 0).

    Newer toolkits also emit architecture-specific suffixes such as ``sm_90a``
    and ``sm_120a``, which are stripped before parsing.
    """

    if not flag.startswith("sm_"):
        return None
    digits = flag[3:].rstrip("abcdefghijklmnopqrstuvwxyz")
    if not digits.isdigit() or len(digits) < 2:
        return None
    return int(digits[:-1]), int(digits[-1])


def execution_path(
    capability: tuple[int, int], arch_list: list[str]
) -> dict[str, Any]:
    """Say in words how this build reaches this device, or that it cannot.

    Three cases have to stay distinguishable in the output file, because a later
    reader will otherwise collapse them.  A device can be served by its own
    native SASS; by SASS for an older minor version of the same major
    architecture, which is binary compatible but not native; or by nothing at
    all, in which case only embedded PTX could rescue it.

    The L40S is the case that misleads.  It is sm_89, and the cu130 build's
    arch list stops at sm_86 within major version 8, so on that runtime the
    L40S runs sm_86 binaries.  That must not be read later as native sm_89.
    """

    major, minor = capability
    device_sm = f"sm_{major}{minor}"
    # The last digit is the minor version and everything before it is the
    # major, so sm_86 is 8.6 while sm_120 is 12.0.  Slicing a fixed number of
    # characters gets Blackwell wrong.
    same_major = sorted(
        (
            parsed[1]
            for parsed in (parse_arch(flag) for flag in arch_list)
            if parsed is not None and parsed[0] == major
        ),
        reverse=True,
    )
    if device_sm in arch_list:
        return {
            "execution_path": "native_sass",
            "execution_note": (
                f"{device_sm} is present in the build's arch list; the device "
                "runs SASS compiled for its own compute capability."
            ),
        }
    if same_major:
        served_by = f"sm_{major}{same_major[0]}"
        return {
            "execution_path": "minor_version_binary_compatibility",
            "execution_served_by": served_by,
            "execution_note": (
                f"{device_sm} is NOT in the build's arch list. The device runs "
                f"{served_by} SASS by CUDA minor-version binary compatibility "
                f"within major architecture {major}. This is compatible but it "
                f"is NOT native {device_sm} execution, and it is not tuned for "
                "this device."
            ),
        }
    return {
        "execution_path": "no_matching_sass",
        "execution_note": (
            f"{device_sm} has no SASS in this build and no other arch of major "
            f"version {major} is present. Execution depends entirely on JIT "
            "from embedded PTX, which requires a compute_* entry in the arch "
            "list."
        ),
    }


# --------------------------------------------------------------------------
# Probe 2: does a real kernel run, and are its numbers right
# --------------------------------------------------------------------------
def probe_kernel(device: torch.device, seed: int) -> dict[str, Any]:
    """Run matmul and conv2d and check them against a float64 CPU reference.

    "It did not crash" is not the bar.  A kernel that silently returns garbage
    on an unsupported architecture would pass a crash-only check, so both
    operations are scored against a reference computed on the host.
    """

    record: dict[str, Any] = {}
    generator = torch.Generator(device="cpu").manual_seed(seed)
    left = torch.randn(512, 384, generator=generator, dtype=torch.float32)
    right = torch.randn(384, 256, generator=generator, dtype=torch.float32)
    reference = (left.to(torch.float64) @ right.to(torch.float64)).to(torch.float32)
    observed = (left.to(device) @ right.to(device)).cpu()
    record["matmul"] = {
        "shape": list(observed.shape),
        "max_abs_error_vs_cpu_float64": float((observed - reference).abs().max()),
        "all_finite": bool(torch.isfinite(observed).all()),
        "digest": tensor_digest(observed),
    }

    image = torch.randn(8, 3, 64, 64, generator=generator, dtype=torch.float32)
    weight = torch.randn(16, 3, 3, 3, generator=generator, dtype=torch.float32)
    conv_reference = nn.functional.conv2d(
        image.to(torch.float64), weight.to(torch.float64), padding=1
    ).to(torch.float32)
    conv_observed = nn.functional.conv2d(
        image.to(device), weight.to(device), padding=1
    ).cpu()
    record["conv2d"] = {
        "shape": list(conv_observed.shape),
        "max_abs_error_vs_cpu_float64": float(
            (conv_observed - conv_reference).abs().max()
        ),
        "all_finite": bool(torch.isfinite(conv_observed).all()),
        "digest": tensor_digest(conv_observed),
    }

    # Reduced precision is where architectures diverge most, and it is what the
    # campaign's models actually run in.  Report it, but do not check on it.
    for name, dtype in (("bfloat16", torch.bfloat16), ("float16", torch.float16)):
        try:
            low = (left.to(device, dtype) @ right.to(device, dtype)).float().cpu()
        except Exception as error:  # noqa: BLE001
            record[f"matmul_{name}"] = f"unavailable: {type(error).__name__}: {error}"
            continue
        record[f"matmul_{name}"] = {
            "max_abs_error_vs_cpu_float64": float((low - reference).abs().max()),
            "all_finite": bool(torch.isfinite(low).all()),
            "digest": tensor_digest(low),
        }

    tolerance = 2e-3
    record["status"] = (
        "pass"
        if (
            record["matmul"]["all_finite"]
            and record["conv2d"]["all_finite"]
            and record["matmul"]["max_abs_error_vs_cpu_float64"] < tolerance
            and record["conv2d"]["max_abs_error_vs_cpu_float64"] < tolerance
        )
        else "fail"
    )
    record["float32_tolerance"] = tolerance
    return record


# --------------------------------------------------------------------------
# Probe 3: strict, weights-only load of an included checkpoint
# --------------------------------------------------------------------------
def module_matching(state: dict[str, torch.Tensor]) -> nn.Module:
    """Build a module whose parameter tree mirrors the checkpoint's own keys.

    ``strict=True`` is only meaningful against a module that claims exactly the
    checkpoint's keys, so the module is derived from the checkpoint rather than
    hand-written.
    """

    root = nn.Module()
    for key, value in state.items():
        *parents, leaf = key.split(".")
        node = root
        for part in parents:
            child = getattr(node, part, None)
            if child is None:
                child = nn.Module()
                node.add_module(part, child)
            node = child
        empty = torch.empty_like(value, device="cpu")
        if value.is_floating_point():
            node.register_parameter(leaf, nn.Parameter(empty))
        else:
            node.register_buffer(leaf, empty)
    return root


def probe_checkpoint_load(
    checkpoint: Path, expected_sha256: str | None
) -> dict[str, Any]:
    record: dict[str, Any] = {"checkpoint": str(checkpoint)}
    observed_sha256 = sha256_file(checkpoint)
    record["checkpoint_sha256"] = observed_sha256
    if expected_sha256 is not None and observed_sha256 != expected_sha256:
        raise ProbeError(
            f"checkpoint digest differs: {observed_sha256} != {expected_sha256}"
        )
    record["checkpoint_sha256_verified"] = expected_sha256 is not None

    # weights_only=True is the whole point: the included checkpoint must restore
    # under the restricted unpickler, with no unrestricted torch.load fallback.
    state = torch.load(checkpoint, map_location="cpu", weights_only=True)
    if not isinstance(state, dict):
        raise ProbeError(f"checkpoint is not a state dict: {type(state).__name__}")
    tensors = {
        key: value for key, value in state.items() if isinstance(value, torch.Tensor)
    }
    if len(tensors) != len(state):
        raise ProbeError("checkpoint holds non-tensor entries")

    module = module_matching(tensors)
    result = module.load_state_dict(tensors, strict=True)
    record["missing_keys"] = list(result.missing_keys)
    record["unexpected_keys"] = list(result.unexpected_keys)
    record["tensor_count"] = len(tensors)
    record["parameter_count"] = int(sum(value.numel() for value in tensors.values()))
    record["dtypes"] = sorted({str(value.dtype) for value in tensors.values()})
    record["state_digest"] = sha256(
        json.dumps(
            {key: tensor_digest(value) for key, value in sorted(tensors.items())},
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()
    record["status"] = (
        "pass" if not result.missing_keys and not result.unexpected_keys else "fail"
    )
    return record


# --------------------------------------------------------------------------
# Probes 4 and 5: a training step that moves, and a resume that matches
# --------------------------------------------------------------------------
def synthetic_task(
    device: torch.device, seed: int
) -> tuple[torch.Tensor, torch.Tensor]:
    generator = torch.Generator(device="cpu").manual_seed(seed)
    features = torch.randn(256, 64, generator=generator)
    truth = torch.randn(64, 8, generator=generator)
    targets = features @ truth + 0.01 * torch.randn(256, 8, generator=generator)
    return features.to(device), targets.to(device)


def build_model(device: torch.device, seed: int) -> nn.Module:
    torch.manual_seed(seed)
    model = nn.Sequential(
        nn.Linear(64, 128),
        nn.ReLU(),
        nn.Linear(128, 64),
        nn.ReLU(),
        nn.Linear(64, 8),
    )
    return model.to(device)


def run_training(
    device: torch.device,
    seed: int,
    steps: int,
    *,
    resume_from: dict[str, Any] | None = None,
) -> tuple[list[float], nn.Module, torch.optim.Optimizer]:
    features, targets = synthetic_task(device, seed)
    model = build_model(device, seed)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    if resume_from is not None:
        model.load_state_dict(resume_from["model"])
        optimizer.load_state_dict(resume_from["optimizer"])
    losses: list[float] = []
    for _ in range(steps):
        optimizer.zero_grad(set_to_none=True)
        loss = nn.functional.mse_loss(model(features), targets)
        loss.backward()
        optimizer.step()
        losses.append(float(loss.detach()))
    return losses, model, optimizer


def probe_training_step(device: torch.device, seed: int) -> dict[str, Any]:
    losses, model, _ = run_training(device, seed, TRAIN_STEPS)
    finite = all(math.isfinite(value) for value in losses)
    record = {
        "steps": TRAIN_STEPS,
        "loss_first": losses[0],
        "loss_last": losses[-1],
        "loss_trajectory": losses,
        "all_finite": finite,
        "decreased": losses[-1] < losses[0],
        "final_parameter_digest": sha256(
            json.dumps(
                {
                    name: tensor_digest(value)
                    for name, value in sorted(model.state_dict().items())
                },
                sort_keys=True,
            ).encode("utf-8")
        ).hexdigest(),
    }
    record["status"] = "pass" if finite and record["decreased"] else "fail"
    return record


def probe_resume(device: torch.device, seed: int, stage: Path) -> dict[str, Any]:
    half = TRAIN_STEPS // 2
    uninterrupted, reference_model, _ = run_training(device, seed, TRAIN_STEPS)

    first_half, model, optimizer = run_training(device, seed, half)
    checkpoint = stage / "resume_state.pt"
    torch.save(
        {
            "model": {k: v.cpu() for k, v in model.state_dict().items()},
            "optimizer": optimizer.state_dict(),
        },
        checkpoint,
    )
    restored = torch.load(checkpoint, map_location=device, weights_only=True)
    second_half, resumed_model, _ = run_training(
        device, seed, half, resume_from=restored
    )

    reference_state = reference_model.state_dict()
    resumed_state = resumed_model.state_dict()
    deltas = [
        float((resumed_state[name].to(torch.float64) - value.to(torch.float64)).abs().max())
        for name, value in reference_state.items()
    ]
    record = {
        "steps_total": TRAIN_STEPS,
        "steps_before_save": half,
        "resume_checkpoint_sha256": sha256_file(checkpoint),
        "uninterrupted_final_loss": uninterrupted[-1],
        "resumed_final_loss": second_half[-1],
        "max_abs_parameter_difference": max(deltas) if deltas else 0.0,
        "loss_absolute_difference": abs(uninterrupted[-1] - second_half[-1]),
        "first_half_matches": first_half == uninterrupted[:half],
    }
    # Same device, same library, same seed: a correct resume is exact.
    record["status"] = (
        "pass"
        if record["max_abs_parameter_difference"] == 0.0
        and record["first_half_matches"]
        else "fail"
    )
    return record


# --------------------------------------------------------------------------
# Probe 6: the numeric fixture the cells are compared on
# --------------------------------------------------------------------------
def probe_numeric_fixture(
    device: torch.device, seed: int, stage: Path
) -> dict[str, Any]:
    """Emit the seeded fixture whose cross-cell difference decides inclusion.

    Raw tensors are written alongside the digests so the comparison between two
    cells is an exact elementwise difference, not a hash equality test that can
    only say "same" or "different".
    """

    generator = torch.Generator(device="cpu").manual_seed(seed)
    outputs: dict[str, torch.Tensor] = {}

    left = torch.randn(256, 256, generator=generator)
    right = torch.randn(256, 256, generator=generator)
    outputs["matmul_fp32"] = (left.to(device) @ right.to(device)).cpu()

    image = torch.randn(4, 8, 32, 32, generator=generator)
    weight = torch.randn(8, 8, 3, 3, generator=generator)
    outputs["conv2d_fp32"] = nn.functional.conv2d(
        image.to(device), weight.to(device), padding=1
    ).cpu()

    outputs["softmax_fp32"] = torch.softmax(
        (left.to(device) * 3.0), dim=-1
    ).cpu()
    outputs["layernorm_fp32"] = nn.functional.layer_norm(
        left.to(device), (256,)
    ).cpu()

    losses, model, _ = run_training(device, seed, TRAIN_STEPS)
    for name, value in model.state_dict().items():
        outputs[f"trained_{name}"] = value.detach().cpu()

    torch.save(outputs, stage / "numeric_fixture.pt")

    # Every input above is drawn on CPU from an explicit seeded generator and
    # only then moved to the device, so the three cells provably start from
    # identical bits and any difference in the outputs is device arithmetic
    # rather than a divergent RNG stream.  These digests are what let a reader
    # confirm that rather than take it on trust.
    initial = build_model(torch.device("cpu"), seed).state_dict()
    inputs = {"features": synthetic_task(torch.device("cpu"), seed)[0]}
    return {
        "seed": seed,
        "initial_parameter_digests": {
            name: tensor_digest(value) for name, value in sorted(initial.items())
        },
        "input_digests": {
            name: tensor_digest(value) for name, value in sorted(inputs.items())
        },
        "tensors": {
            name: {
                "shape": list(value.shape),
                "dtype": str(value.dtype),
                "digest": tensor_digest(value),
                "mean": float(value.to(torch.float64).mean()),
                "absmax": float(value.to(torch.float64).abs().max()),
            }
            for name, value in sorted(outputs.items())
        },
        "training_loss_trajectory": losses,
        "fixture_path": "numeric_fixture.pt",
        "fixture_sha256": sha256_file(stage / "numeric_fixture.pt"),
        "status": "pass",
    }


# --------------------------------------------------------------------------
# Cross-cell comparison
# --------------------------------------------------------------------------
def compare_fixtures(left: Path, right: Path) -> dict[str, Any]:
    """Elementwise difference between two cells' probe-6 fixtures."""

    left_tensors = torch.load(left, map_location="cpu", weights_only=True)
    right_tensors = torch.load(right, map_location="cpu", weights_only=True)
    shared = sorted(set(left_tensors) & set(right_tensors))
    per_tensor: dict[str, Any] = {}
    worst_absolute = 0.0
    worst_relative = 0.0
    for name in shared:
        a = left_tensors[name].to(torch.float64)
        b = right_tensors[name].to(torch.float64)
        if a.shape != b.shape:
            per_tensor[name] = {"status": "shape_differs"}
            continue
        absolute = (a - b).abs()
        scale = torch.maximum(a.abs(), b.abs()).clamp_min(1e-12)
        relative = absolute / scale
        entry = {
            "max_abs_difference": float(absolute.max()),
            "max_rel_difference": float(relative.max()),
            "mean_abs_difference": float(absolute.mean()),
            "bitwise_identical": bool(torch.equal(a, b)),
        }
        per_tensor[name] = entry
        worst_absolute = max(worst_absolute, entry["max_abs_difference"])
        worst_relative = max(worst_relative, entry["max_rel_difference"])
    return {
        "schema_version": COMPARE_SCHEMA,
        "left_fixture": str(left),
        "right_fixture": str(right),
        "left_fixture_sha256": sha256_file(left),
        "right_fixture_sha256": sha256_file(right),
        "tensors_compared": len(shared),
        "left_only": sorted(set(left_tensors) - set(right_tensors)),
        "right_only": sorted(set(right_tensors) - set(left_tensors)),
        "max_abs_difference": worst_absolute,
        "max_rel_difference": worst_relative,
        "all_bitwise_identical": all(
            entry.get("bitwise_identical") is True for entry in per_tensor.values()
        ),
        "per_tensor": per_tensor,
    }


# --------------------------------------------------------------------------
def guarded(action: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    """Run one probe, converting any failure into a recorded result.

    A build with no kernel for the device is an expected terminal disposition,
    not a crash to propagate.  The campaign needs that answer written down as
    an output file, so no probe is allowed to abort the run.
    """

    try:
        return action()
    except Exception as error:  # noqa: BLE001 - the failure is the measurement
        return {
            "status": "error",
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(limit=8),
        }


def run_all(
    cell_id: str, checkpoint: Path, checkpoint_sha256: str | None, seed: int, stage: Path
) -> dict[str, Any]:
    determinism = pin_determinism(seed)
    probes: dict[str, Any] = {}
    probes["probe_1_import_architecture"] = guarded(probe_import_architecture)

    skip = None
    if probes["probe_1_import_architecture"].get("cuda_available") is not True:
        skip = "skipped_no_visible_cuda_device"

    device = torch.device("cuda:0")
    if skip is None:
        # Probe 2 first: if the build has no kernel for this architecture, every
        # later probe would fail for the same single reason, and saying so once
        # is more useful than saying it five times.
        probes["probe_2_kernel"] = guarded(lambda: probe_kernel(device, seed))
        if probes["probe_2_kernel"]["status"] != "pass":
            skip = "skipped_no_working_kernel"
    else:
        probes["probe_2_kernel"] = {"status": skip}

    remaining: list[tuple[str, Callable[[], dict[str, Any]]]] = [
        (
            "probe_3_checkpoint_load",
            lambda: probe_checkpoint_load(checkpoint, checkpoint_sha256),
        ),
        ("probe_4_training_step", lambda: probe_training_step(device, seed)),
        ("probe_5_resume", lambda: probe_resume(device, seed, stage)),
        ("probe_6_numeric_fixture", lambda: probe_numeric_fixture(device, seed, stage)),
    ]
    for name, action in remaining:
        probes[name] = {"status": skip} if skip is not None else guarded(action)

    statuses = {name: value.get("status") for name, value in probes.items()}
    return {
        "schema_version": PROBE_SCHEMA,
        "cell_id": cell_id,
        "seed": seed,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", "interactive"),
        "node": platform.node(),
        "determinism_knobs": determinism,
        "outcomes_read": False,
        "project_data_read": False,
        "metrics_calculated": False,
        "probe_statuses": statuses,
        "status": "pass" if set(statuses.values()) == {"pass"} else "fail",
        "probes": probes,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-id")
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--checkpoint-sha256")
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare", nargs=2, type=Path, metavar=("LEFT", "RIGHT"))
    args = parser.parse_args()

    if args.compare is not None:
        receipt = compare_fixtures(args.compare[0], args.compare[1])
    else:
        for name in ("cell_id", "checkpoint", "output"):
            if getattr(args, name) is None:
                raise SystemExit(f"--{name.replace('_', '-')} is required")
        args.output.mkdir(parents=True, exist_ok=True)
        receipt = run_all(
            args.cell_id,
            args.checkpoint.resolve(strict=True),
            args.checkpoint_sha256,
            args.seed,
            args.output,
        )

    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        name = "comparison.json" if args.compare is not None else "probe.json"
        (args.output / name).write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if receipt.get("status", "pass") == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
