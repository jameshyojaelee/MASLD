#!/usr/bin/env python3
"""Strictly restore Scooby and probe its checkpoint decoder on synthetic inputs."""

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


SCHEMA = "masld-bench-scooby-epicardioids-checkpoint-decoder-probe-v1"


class ScoobyCheckpointDecoderProbeError(RuntimeError):
    """Raised when strict restore or the synthetic checkpoint probe differs."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ScoobyCheckpointDecoderProbeError(message)


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
        raise ScoobyCheckpointDecoderProbeError(f"invalid {label}: {error}") from error
    require(isinstance(value, dict), f"{label} is not an object")
    return value


def verify_artifact(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    path = root / str(binding["path"])
    try:
        verify_frozen_tree(path)
    except ArtifactError as error:
        raise ScoobyCheckpointDecoderProbeError(f"{label} artifact differs: {error}") from error
    require(digest(path / "ARTIFACTS.json") == binding["artifacts_sha256"], f"{label} artifact identity differs")
    return path


def validate_task_boundary(config: Mapping[str, Any]) -> None:
    task = config["task_scope"]
    policy = config["probe_policy"]
    require(task["observed_atac_in_context"], "observed-ATAC provenance was lost")
    require(not task["rna_conditioned_atac_eligible"], "checkpoint cannot enter RNA-conditioned ATAC lane")
    require(not task["liver_query_encoder_available"], "unreleased liver query encoder was asserted")
    require(not task["sealed_inference_allowed"], "probe cannot authorize sealed inference")
    require(not task["universal_or_champion_claim_allowed"], "probe cannot authorize a champion claim")
    require(policy["safetensors_only"], "probe must use safetensors only")
    require(policy["strict_restore_required"], "strict state restore is required")
    require(policy["synthetic_decoder_forward_allowed"], "synthetic decoder probe was not authorized")
    require(not policy["full_sequence_forward_allowed"], "full sequence forward is outside this gate")
    require(not policy["native_cell_context_claim_allowed"], "synthetic context is not a native cell context")
    require(not policy["project_data_read_allowed"] and not policy["sealed_data_read_allowed"], "data firewall differs")


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


def runtime_probe(root: Path, config: Mapping[str, Any], source_dir: Path, checkpoint: Path) -> dict[str, Any]:
    import torch
    from borzoi_pytorch.config_borzoi import BorzoiConfig
    from safetensors.torch import load_file

    expected_runtime = config["runtime"]
    versions = {
        "torch": torch.__version__,
        "borzoi_pytorch": importlib.metadata.version("borzoi-pytorch"),
        "transformers": importlib.metadata.version("transformers"),
        "safetensors": importlib.metadata.version("safetensors"),
        "einops": importlib.metadata.version("einops"),
    }
    require(versions == {key: expected_runtime[key] for key in versions}, "runtime versions differ")
    require(checkpoint.suffix == ".safetensors", "checkpoint is not safetensors")
    require(checkpoint.stat().st_size == config["checkpoint"]["size_bytes"], "checkpoint size differs")
    require(digest(checkpoint) == config["checkpoint"]["sha256"], "checkpoint hash differs before restore")
    native_config = load_json(root / config["authorities"]["native_config"]["path"], label="native config")
    constructor = config["constructor"]
    scooby = load_scooby(source_dir / "scooby/modeling/scooby.py")
    torch.manual_seed(config["fixture"]["seed"])
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
    model.eval()

    cell_context = torch.linspace(-1.0, 1.0, steps=constructor["cell_emb_dim"], dtype=torch.float32).reshape(1, 1, -1)
    sequence_embedding = torch.linspace(
        -0.25,
        0.25,
        steps=constructor["embedding_dim"] * config["fixture"]["synthetic_bins"],
        dtype=torch.float32,
    ).reshape(1, constructor["embedding_dim"], config["fixture"]["synthetic_bins"])
    with torch.inference_mode():
        weights_a, biases_a = model.forward_cell_embs_only(cell_context)
        profile_a = model.forward_convs_on_emb(sequence_embedding, weights_a, biases_a)
        weights_b, biases_b = model.forward_cell_embs_only(cell_context)
        profile_b = model.forward_convs_on_emb(sequence_embedding, weights_b, biases_b)
    require(torch.equal(weights_a, weights_b) and torch.equal(biases_a, biases_b), "checkpoint decoder repeat differs")
    require(torch.equal(profile_a, profile_b), "checkpoint profile repeat differs")
    require(tuple(weights_a.shape) == (1, 3, 1920, 1), "checkpoint decoder weight shape differs")
    require(tuple(biases_a.shape) == (1, 3), "checkpoint decoder bias shape differs")
    require(tuple(profile_a.shape) == (1, config["fixture"]["synthetic_bins"], 3), "checkpoint profile shape differs")
    require(torch.isfinite(profile_a).all().item() and torch.all(profile_a > 0).item(), "checkpoint profile numeric contract differs")
    return {
        "runtime_versions": versions,
        "safe_loader": "safetensors.torch.load_file",
        "strict_restore": True,
        "missing_keys": [],
        "unexpected_keys": [],
        "state_dict_tensors": constructor["expected_state_dict_tensors"],
        "state_dict_elements_including_buffers": constructor["expected_state_dict_elements_including_buffers"],
        "synthetic_context": tensor_receipt(cell_context),
        "synthetic_sequence_embedding": tensor_receipt(sequence_embedding),
        "decoder_weights": tensor_receipt(weights_a),
        "decoder_biases": tensor_receipt(biases_a),
        "profile": tensor_receipt(profile_a),
        "deterministic_repeat_bit_identical": True,
        "checkpoint_deserialized_with_safetensors": True,
        "checkpoint_decoder_forward_executed": True,
        "full_sequence_forward_executed": False,
        "native_cell_context_used": False,
        "project_data_read": False,
        "sealed_data_read": False,
    }


def audit(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    config = load_json(config_path.resolve(strict=True), label="decoder probe config")
    require(config.get("schema_version") == SCHEMA, "decoder probe config identity differs")
    require(not output.exists() and not output.is_symlink(), "refusing to overwrite probe output")
    validate_task_boundary(config)
    acquisition = verify_artifact(root, config["authorities"]["checkpoint_artifact"], label="checkpoint")
    schema = verify_artifact(root, config["authorities"]["runtime_schema_artifact"], label="runtime schema")
    acquisition_receipt = load_json(acquisition / "audit/receipt.json", label="acquisition receipt")
    require(acquisition_receipt["status"] == "pass_exact_checkpoint_acquired_safe_restore_allowed", "acquisition did not authorize restore")
    require(acquisition_receipt["checkpoint"]["sha256"] == config["checkpoint"]["sha256"], "acquired checkpoint identity differs")
    for name in ("native_config", "official_source"):
        binding = config["authorities"][name]
        require(digest(root / binding["path"]) == binding["sha256"], f"{name} authority differs")
    probe = runtime_probe(root, config, schema / "sources", acquisition / "checkpoint/model.safetensors")
    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": "masld-bench-scooby-epicardioids-checkpoint-decoder-probe-receipt-v1",
        "status": "pass_checkpoint_decoder_executable_full_sequence_and_native_context_still_blocked",
        "config_sha256": digest(config_path),
        "checkpoint": config["checkpoint"],
        "license": config["license"],
        "task_scope": config["task_scope"],
        "runtime_probe": probe,
        "checkpoint_payload_downloaded": True,
        "checkpoint_deserialized_with_safetensors": True,
        "checkpoint_decoder_forward_executed": True,
        "full_sequence_forward_executed": False,
        "native_cell_context_used": False,
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
