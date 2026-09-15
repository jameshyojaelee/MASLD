#!/usr/bin/env python3
"""Run a bounded full-backbone Scooby checkpoint forward on synthetic inputs."""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree, write_json_exclusive
from scripts.audit_scooby_epicardioids_runtime_schema_gate import load_scooby


SCHEMA = "masld-bench-scooby-epicardioids-full-backbone-synthetic-probe-v1"


class ScoobyFullBackboneProbeError(RuntimeError):
    """Raised when the full-backbone synthetic inclusion differs."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ScoobyFullBackboneProbeError(message)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    require(path.is_file() and not path.is_symlink(), f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScoobyFullBackboneProbeError(f"invalid {label}: {error}") from error
    require(isinstance(value, dict), f"{label} is not an object")
    return value


def verify_artifact(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    path = root / str(binding["path"])
    try:
        verify_frozen_tree(path)
    except ArtifactError as error:
        raise ScoobyFullBackboneProbeError(f"{label} artifact differs: {error}") from error
    require(digest(path / "ARTIFACTS.json") == binding["artifacts_sha256"], f"{label} identity differs")
    return path


def validate_boundary(config: Mapping[str, Any]) -> None:
    task = config["task_scope"]
    policy = config["probe_policy"]
    require(task["synthetic_architecture_admission_only"], "probe is not synthetic-only")
    require(not task["native_cell_context_used"], "synthetic context cannot be called native")
    require(not task["biological_task_executed"], "probe cannot execute a biological task")
    require(task["observed_atac_in_released_context"], "released context provenance was lost")
    require(not task["rna_conditioned_atac_eligible"], "observed-ATAC model cannot enter RNA-conditioned lane")
    require(not task["coordinate_compatible_liver_query_supported"], "coordinate-compatible liver query was asserted")
    require(not task["sealed_inference_allowed"] and not task["universal_or_champion_claim_allowed"], "claim firewall differs")
    require(policy["safetensors_only"] and policy["strict_restore_required"], "safe strict restore is required")
    require(policy["full_sequence_synthetic_forward_allowed"], "full synthetic sequence forward was not authorized")
    require(policy["deterministic_repeat_required"], "deterministic repeat is required")
    require(not policy["project_data_read_allowed"] and not policy["sealed_data_read_allowed"], "data firewall differs")
    require(not policy["model_fit_or_adaptation_allowed"], "adaption is outside admission")


def tensor_receipt(tensor: Any) -> dict[str, Any]:
    import numpy as np

    array = tensor.detach().cpu().numpy().astype("<f4", copy=False)
    return {
        "shape": list(array.shape),
        "dtype": "float32_le",
        "sha256": hashlib.sha256(array.tobytes(order="C")).hexdigest(),
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
        "mean": float(np.mean(array, dtype=np.float64)),
        "finite": bool(np.isfinite(array).all()),
    }


def run_probe(root: Path, config: Mapping[str, Any], schema: Path, checkpoint: Path) -> dict[str, Any]:
    import torch
    import torch.nn.functional as torch_f
    from borzoi_pytorch.config_borzoi import BorzoiConfig
    from safetensors.torch import load_file

    require(torch.cuda.is_available(), "CUDA is unavailable")
    expected_runtime = config["runtime"]
    versions = {
        "torch": torch.__version__,
        "borzoi_pytorch": importlib.metadata.version("borzoi-pytorch"),
        "transformers": importlib.metadata.version("transformers"),
        "safetensors": importlib.metadata.version("safetensors"),
        "einops": importlib.metadata.version("einops"),
    }
    require(versions == {key: expected_runtime[key] for key in versions}, "runtime versions differ")
    device = torch.device("cuda:0")
    gpu_name = torch.cuda.get_device_name(device)
    capability = list(torch.cuda.get_device_capability(device))
    require(expected_runtime["required_gpu_name_contains"] in gpu_name, "GPU model differs")
    require(capability == expected_runtime["required_compute_capability"], "GPU compute capability differs")
    torch.manual_seed(config["fixture"]["seed"])
    torch.cuda.manual_seed_all(config["fixture"]["seed"])
    torch.use_deterministic_algorithms(True)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

    for name in ("native_config", "official_source"):
        binding = config["authorities"][name]
        require(digest(root / binding["path"]) == binding["sha256"], f"{name} authority differs")
    native_config = load_json(root / config["authorities"]["native_config"]["path"], label="native config")
    constructor = config["constructor"]
    scooby = load_scooby(schema / "sources/scooby/modeling/scooby.py")
    model = scooby.Scooby(
        BorzoiConfig(**native_config),
        cell_emb_dim=constructor["cell_emb_dim"],
        embedding_dim=constructor["embedding_dim"],
        n_tracks=constructor["n_tracks"],
        disable_cache=constructor["disable_cache"],
        use_transform_borzoi_emb=constructor["use_transform_borzoi_emb"],
    )
    state = load_file(str(checkpoint), device="cpu")
    require(len(state) == constructor["expected_state_dict_tensors"], "loaded state tensor count differs")
    require(sum(tensor.numel() for tensor in state.values()) == constructor["expected_state_dict_elements_including_buffers"], "loaded state element count differs")
    restored = model.load_state_dict(state, strict=True)
    require(not restored.missing_keys and not restored.unexpected_keys, "strict restore returned incompatible keys")
    del state
    gc.collect()
    model.eval().to(device)

    length = config["fixture"]["sequence_length_bp"]
    base_indices = torch.arange(length, dtype=torch.int64).remainder(4)
    sequence_cpu = torch_f.one_hot(base_indices, num_classes=4).to(torch.float32).unsqueeze(0)
    context_cpu = torch.linspace(-1.0, 1.0, steps=constructor["cell_emb_dim"], dtype=torch.float32).reshape(1, 1, -1)
    sequence = sequence_cpu.to(device)
    context = context_cpu.to(device)
    torch.cuda.reset_peak_memory_stats(device)
    with torch.inference_mode():
        output_a = model(sequence, context)
        torch.cuda.synchronize(device)
        output_b = model(sequence, context)
        torch.cuda.synchronize(device)
    require(tuple(output_a.shape) == tuple(config["fixture"]["expected_output_shape"]), "full-backbone output shape differs")
    require(torch.equal(output_a, output_b), "full-backbone deterministic repeat differs")
    require(torch.isfinite(output_a).all().item() and torch.all(output_a > 0).item(), "full-backbone numeric contract differs")
    output_receipt = tensor_receipt(output_a)
    return {
        "runtime_versions": versions,
        "cuda_runtime": torch.version.cuda,
        "gpu_name": gpu_name,
        "gpu_compute_capability": capability,
        "gpu_total_memory_bytes": torch.cuda.get_device_properties(device).total_memory,
        "gpu_peak_allocated_bytes": torch.cuda.max_memory_allocated(device),
        "safe_loader": "safetensors.torch.load_file",
        "strict_restore": True,
        "missing_keys": [],
        "unexpected_keys": [],
        "state_dict_tensors": constructor["expected_state_dict_tensors"],
        "state_dict_elements_including_buffers": constructor["expected_state_dict_elements_including_buffers"],
        "synthetic_sequence": tensor_receipt(sequence_cpu),
        "synthetic_context": tensor_receipt(context_cpu),
        "output": output_receipt,
        "deterministic_repeat_bit_identical": True,
        "checkpoint_deserialized_with_safetensors": True,
        "full_sequence_forward_executed": True,
        "native_cell_context_used": False,
        "biological_task_executed": False,
        "model_fitted_or_adapted": False,
        "project_data_read": False,
        "sealed_data_read": False,
    }


def audit(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    config = load_json(config_path.resolve(strict=True), label="full-backbone config")
    require(config.get("schema_version") == SCHEMA, "full-backbone config identity differs")
    require(not output.exists() and not output.is_symlink(), "refusing to overwrite full-backbone output")
    validate_boundary(config)
    checkpoint_artifact = verify_artifact(root, config["authorities"]["checkpoint_artifact"], label="checkpoint")
    decoder = verify_artifact(root, config["authorities"]["decoder_probe_artifact"], label="decoder probe")
    schema = verify_artifact(root, config["authorities"]["runtime_schema_artifact"], label="runtime schema")
    decoder_receipt = load_json(decoder / "audit/receipt.json", label="decoder probe receipt")
    require(decoder_receipt["status"] == "pass_checkpoint_decoder_executable_full_sequence_and_native_context_still_blocked", "decoder gate did not authorize backbone probe")
    checkpoint = checkpoint_artifact / "checkpoint/model.safetensors"
    require(checkpoint.stat().st_size == config["checkpoint"]["size_bytes"], "checkpoint size differs")
    require(digest(checkpoint) == config["checkpoint"]["sha256"], "checkpoint hash differs")
    probe = run_probe(root, config, schema, checkpoint)
    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": "masld-bench-scooby-epicardioids-full-backbone-synthetic-probe-receipt-v1",
        "status": "pass_full_backbone_synthetic_native_liver_context_still_blocked",
        "config_sha256": digest(config_path),
        "checkpoint": config["checkpoint"],
        "task_scope": config["task_scope"],
        "runtime_probe": probe,
        "full_sequence_forward_executed": True,
        "native_cell_context_used": False,
        "biological_task_executed": False,
        "model_fitted_or_adapted": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "next_gate": config["probe_policy"]["next_if_pass"],
    }
    write_json_exclusive(output / "receipt.json", receipt, mode=0o640)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.root, args.config, args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
