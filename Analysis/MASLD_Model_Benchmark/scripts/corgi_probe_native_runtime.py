#!/usr/bin/env python3
"""Strict Corgi restore and shortened-window GPU compatibility probe.

This is not a native numeric-parity test. It loads only verified tensor state,
uses the exact released architecture dimensions, and shortens only the spatial
input/crop so forward and backward kernels can be checked without a full-window
benchmark run.
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
from pathlib import Path
import random
import sys
import types
from typing import Any

import torch


REGULAR_SHA256 = "cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4"
PLUS_SHA256 = "01855614bfeffa72d5bcaaa28729c022ccfbc8f5fdcc0add8e04ba4533e94ccd"


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--regular", type=Path, required=True)
    parser.add_argument("--plus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def _load_released_model_source(source: Path) -> tuple[type, type, dict[str, Any]]:
    package_dir = source / "corgi"
    package = types.ModuleType("corgi")
    package.__path__ = [str(package_dir)]
    package.__package__ = "corgi"
    sys.modules["corgi"] = package
    loaded: dict[str, Any] = {}
    for name in ("modules", "model", "config"):
        path = package_dir / f"{name}.py"
        spec = importlib.util.spec_from_file_location(f"corgi.{name}", path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"cannot construct source spec: {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        loaded[name] = module
    config = dict(loaded["config"].config_corgi)
    return loaded["model"].Corgi, loaded["model"].CorgiPlus, config


def _load_state(path: Path) -> dict[str, torch.Tensor]:
    checkpoint = torch.load(str(path), map_location="cpu", weights_only=True, mmap=True)
    if not isinstance(checkpoint, dict) or "model_state_dict" not in checkpoint:
        raise RuntimeError(
            f"checkpoint wrapper differs: {sorted(checkpoint) if isinstance(checkpoint, dict) else type(checkpoint)}"
        )
    state = checkpoint["model_state_dict"]
    if not isinstance(state, dict) or not state:
        raise RuntimeError("checkpoint state_dict is empty or not a mapping")
    if not all(isinstance(key, str) and isinstance(value, torch.Tensor) for key, value in state.items()):
        raise RuntimeError("checkpoint state_dict contains non-tensor records")
    return state


def _network_namespace_record() -> dict[str, Any]:
    route_path = Path("/proc/net/route")
    lines = route_path.read_text(encoding="utf-8").splitlines()
    routes = [line.split() for line in lines[1:] if line.strip()]
    non_loopback = [fields for fields in routes if fields and fields[0] != "lo"]
    if non_loopback:
        raise RuntimeError("runtime probe has a non-loopback network route")
    return {
        "network_namespace": "none",
        "non_loopback_route_count": 0,
        "singularity_no_network_asserted": True,
    }


def _base_config(config: dict[str, Any]) -> dict[str, Any]:
    exact = dict(config)
    exact["output_central_bins"] = 32
    if exact["dim"] != 1536 or exact["input_trans_regulators"] != 2891:
        raise RuntimeError("released architecture dimensions differ")
    if exact["output_channels"] != 22 or exact["heads"] != 8 or exact["gqa"] != 2:
        raise RuntimeError("released output/attention contract differs")
    return exact


def _one_hot(batch: int, length: int, device: torch.device) -> torch.Tensor:
    bases = torch.arange(length, device=device).remainder(4)
    bases = bases.unsqueeze(0).expand(batch, -1)
    return torch.nn.functional.one_hot(bases, num_classes=4).to(torch.bfloat16)


def _reverse_complement(sequence: torch.Tensor) -> torch.Tensor:
    return torch.flip(sequence, dims=(1,)).index_select(
        2, torch.tensor([3, 2, 1, 0], device=sequence.device)
    )


def _probe_one(
    *,
    model_class: type,
    config: dict[str, Any],
    checkpoint_path: Path,
    expected_sha256: str,
    expected_keys: int,
    expected_numel: int,
    plus: bool,
    device: torch.device,
) -> dict[str, Any]:
    state = _load_state(checkpoint_path)
    if len(state) != expected_keys:
        raise RuntimeError(f"state key count differs: {len(state)}")
    total_numel = sum(int(tensor.numel()) for tensor in state.values())
    if total_numel != expected_numel:
        raise RuntimeError(f"state value count differs: {total_numel}")
    if plus:
        config["corgiplus_aux_input_dim"] = 26
        aux_shape = list(state["aux_encoder.0.weight"].shape)
        if aux_shape != [128, 26, 1]:
            raise RuntimeError(f"Corgi+ auxiliary state differs: {aux_shape}")
    direct_output_head = "output_head.weight" in state
    sequential_output_head = "output_head.0.weight" in state
    if direct_output_head == sequential_output_head:
        raise RuntimeError("checkpoint output-head key layout is ambiguous")
    config["final_softplus"] = sequential_output_head
    model = model_class(config)
    incompatible = model.load_state_dict(state, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise RuntimeError(f"strict restore differs: {incompatible}")
    del state
    gc.collect()

    model = model.to(device=device)
    model.eval()
    sequence = _one_hot(1, 4096, device)
    context = torch.linspace(-1.0, 1.0, 2891, device=device, dtype=torch.bfloat16).unsqueeze(0)
    inputs: tuple[torch.Tensor, ...]
    if plus:
        aux = torch.linspace(
            -0.5, 0.5, 64 * 26, device=device, dtype=torch.bfloat16
        ).reshape(1, 64, 26)
        inputs = (sequence, aux, context)
    else:
        inputs = (sequence, context)

    model.crop.target_length = 6144
    torch.cuda.reset_peak_memory_stats(device)
    full_sequence = _one_hot(1, 524288, device)
    if plus:
        full_aux = torch.zeros((1, 8192, 26), device=device, dtype=torch.bfloat16)
        full_inputs = (full_sequence, full_aux, context)
    else:
        full_inputs = (full_sequence, context)
    with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        full_output = model(*full_inputs)
    if list(full_output.shape) != [1, 22, 6144] or not torch.isfinite(full_output).all():
        raise RuntimeError(f"full-window output contract differs: {list(full_output.shape)}")
    full_peak_memory_bytes = int(torch.cuda.max_memory_allocated(device))
    del full_output, full_inputs, full_sequence
    if plus:
        del full_aux
    torch.cuda.empty_cache()
    model.crop.target_length = 32

    optimizer = torch.optim.SGD(model.parameters(), lr=1.0e-8)
    optimizer.zero_grad(set_to_none=True)
    with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        output = model(*inputs)
        loss = output.float().square().mean()
    if list(output.shape) != [1, 22, 32] or not torch.isfinite(output).all():
        raise RuntimeError(f"shortened-window output contract differs: {list(output.shape)}")
    loss.backward()
    gradients = {
        "conv_0": model.conv_0.weight.grad,
        "output_head": (
            model.output_head[0].weight.grad
            if sequential_output_head
            else model.output_head.weight.grad
        ),
    }
    if not all(value is not None and torch.isfinite(value).all() for value in gradients.values()):
        raise RuntimeError("required finite gradients were not produced")
    optimizer.step()
    optimizer.zero_grad(set_to_none=True)

    with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        repeat_a = model(*inputs)
        repeat_b = model(*inputs)
        rc_inputs = (_reverse_complement(sequence),) + inputs[1:]
        reverse_complement = model(*rc_inputs)
    repeat_max_abs = float((repeat_a.float() - repeat_b.float()).abs().max().item())
    if repeat_max_abs > 1.0e-3:
        raise RuntimeError(f"repeat inference tolerance failed: {repeat_max_abs}")
    if list(reverse_complement.shape) != [1, 22, 32] or not torch.isfinite(reverse_complement).all():
        raise RuntimeError("reverse-complement shortened-window shape/finite check failed")

    result = {
        "checkpoint_sha256_expected": expected_sha256,
        "state_key_count": expected_keys,
        "state_numel": expected_numel,
        "strict_restore": True,
        "output_head_layout": (
            "conv1d_then_softplus" if sequential_output_head else "direct_conv1d"
        ),
        "final_softplus_checkpoint_derived": sequential_output_head,
        "missing_keys": [],
        "unexpected_keys": [],
        "probe_input_bp": 4096,
        "probe_output_bins": 32,
        "probe_output_shape": [1, 22, 32],
        "probe_dtype": "bfloat16_autocast_with_float32_master_weights",
        "finite_forward": True,
        "finite_backward": True,
        "optimizer_step": "SGD_lr_1e-8",
        "repeat_max_abs": repeat_max_abs,
        "repeat_tolerance": 1.0e-3,
        "reverse_complement_shape_and_finite_only": True,
        "full_524288_bp_forward_performed": True,
        "full_524288_bp_output_shape": [1, 22, 6144],
        "full_524288_bp_output_finite": True,
        "full_524288_bp_peak_memory_bytes": full_peak_memory_bytes,
        "native_numeric_parity_performed": False,
    }
    if plus:
        result["auxiliary_input_width"] = 26
        result["released_input_bundle_complete"] = False

    del optimizer, model, sequence, context, inputs, output, loss
    gc.collect()
    torch.cuda.empty_cache()
    return result


def main() -> None:
    args = _arguments()
    if not torch.cuda.is_available():
        raise SystemExit("CUDA is unavailable")
    if torch.cuda.device_count() != 1:
        raise SystemExit(f"expected exactly one visible GPU, found {torch.cuda.device_count()}")
    seed = 20260824
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    device = torch.device("cuda:0")
    Corgi, CorgiPlus, released_config = _load_released_model_source(args.source)
    config = _base_config(released_config)
    result = {
        "schema_version": "masld-bench-corgi-native-runtime-probe-v1",
        "status": "pass",
        "seed": seed,
        "torch": torch.__version__,
        "torch_cuda": torch.version.cuda,
        "gpu_name": torch.cuda.get_device_name(device),
        "gpu_capability": list(torch.cuda.get_device_capability(device)),
        "source_loaded_without_package_init": True,
        "released_high_level_loader_used": False,
        "safe_loader": "torch.load(weights_only=True,map_location=cpu,mmap=True)",
        "outcomes_exposed": False,
        "sealed_features_exposed": False,
        "network": _network_namespace_record(),
    }
    result["regular"] = _probe_one(
        model_class=Corgi,
        config=dict(config),
        checkpoint_path=args.regular,
        expected_sha256=REGULAR_SHA256,
        expected_keys=177,
        expected_numel=195870252,
        plus=False,
        device=device,
    )
    result["plus"] = _probe_one(
        model_class=CorgiPlus,
        config=dict(config),
        checkpoint_path=args.plus,
        expected_sha256=PLUS_SHA256,
        expected_keys=190,
        expected_numel=198664365,
        plus=True,
        device=device,
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
