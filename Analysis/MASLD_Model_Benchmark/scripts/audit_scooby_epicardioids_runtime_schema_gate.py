#!/usr/bin/env python3
"""Validate Scooby Epicardioids runtime schema before downloading its checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import sys
import types
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-scooby-epicardioids-runtime-schema-gate-v1"


class ScoobySchemaGateError(RuntimeError):
    """Raised when the checkpoint schema and executable constructor differ."""


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ScoobySchemaGateError(message)


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    require(path.is_file() and not path.is_symlink(), f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScoobySchemaGateError(f"invalid {label}: {error}") from error
    require(isinstance(value, dict), f"{label} is not an object")
    return value


def validate_authorities(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    authorities = config["frozen_authorities"]
    artifact = root / authorities["bounded_header_artifact"]["path"]
    try:
        verify_frozen_tree(artifact)
    except ArtifactError as error:
        raise ScoobySchemaGateError(f"bounded header artifact differs: {error}") from error
    require(
        digest(artifact / "ARTIFACTS.json") == authorities["bounded_header_artifact"]["artifacts_sha256"],
        "bounded header artifact identity differs",
    )
    result: dict[str, Any] = {}
    for name in ("checkpoint_contract", "header", "config", "model_card"):
        binding = authorities[name]
        path = root / binding["path"]
        require(digest(path) == binding["sha256"], f"frozen authority differs: {name}")
        result[name] = {"path": binding["path"], "sha256": binding["sha256"]}
    checkpoint_contract = load_json(root / authorities["checkpoint_contract"]["path"], label="checkpoint contract")
    checkpoint = config["official_checkpoint"]
    require(
        checkpoint_contract["checkpoint"]["revision"] == checkpoint["revision"]
        and checkpoint_contract["checkpoint"]["sha256"] == checkpoint["sha256"]
        and checkpoint_contract["checkpoint"]["size_bytes"] == checkpoint["size_bytes"]
        and checkpoint_contract["weight_license"].startswith("MIT_DECLARED_WITH_UNRESOLVED"),
        "checkpoint contract differs",
    )
    return result


def validate_sources(source_dir: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    records = {}
    for relative, expected in config["official_code"]["files"].items():
        path = source_dir / relative
        require(path.is_file() and not path.is_symlink(), f"source absent: {relative}")
        require(path.stat().st_size == expected["size_bytes"], f"source size differs: {relative}")
        require(digest(path) == expected["sha256"], f"source hash differs: {relative}")
        records[relative] = expected
    hf = load_json(source_dir / "hf_api.json", label="Hugging Face metadata")
    checkpoint = config["official_checkpoint"]
    require(
        hf.get("sha") == checkpoint["revision"]
        and hf.get("lastModified") == checkpoint["last_modified"]
        and hf.get("gated") is False
        and hf.get("private") is False
        and (hf.get("cardData") or {}).get("license") == checkpoint["card_license"],
        "live checkpoint identity or license differs",
    )
    parameters = ((hf.get("safetensors") or {}).get("parameters") or {}).get("F32")
    require(parameters == checkpoint["declared_parameter_count"], "live parameter count differs")
    return {"files": records, "hf_revision": hf["sha"], "hf_card_license": hf["cardData"]["license"]}


def install_peft_inference_placeholder() -> None:
    require("peft" not in sys.modules, "peft was imported before the inference shim")
    module = types.ModuleType("peft")

    class BlockedTrainingApi:
        def __init__(self, *args: object, **kwargs: object):
            raise ScoobySchemaGateError("PEFT training APIs are blocked in the inference-only gate")

    def blocked_get_peft_model(*args: object, **kwargs: object) -> object:
        raise ScoobySchemaGateError("PEFT mutation is blocked in the inference-only gate")

    module.LoraConfig = BlockedTrainingApi
    module.get_peft_model = blocked_get_peft_model
    sys.modules["peft"] = module


def load_scooby(source: Path) -> Any:
    install_peft_inference_placeholder()
    spec = importlib.util.spec_from_file_location("scooby_epicardioids_official", source)
    require(spec is not None and spec.loader is not None, "official Scooby source cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def runtime_probe(root: Path, source_dir: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    import numpy as np
    import torch
    from borzoi_pytorch.config_borzoi import BorzoiConfig

    expected_runtime = config["runtime"]
    versions = {
        "torch": torch.__version__,
        "borzoi_pytorch": importlib.metadata.version("borzoi-pytorch"),
        "transformers": importlib.metadata.version("transformers"),
        "safetensors": importlib.metadata.version("safetensors"),
        "einops": importlib.metadata.version("einops"),
    }
    require(versions == {key: expected_runtime[key] for key in versions}, "runtime versions differ")
    scooby = load_scooby(source_dir / "scooby/modeling/scooby.py")
    native_config = load_json(root / config["frozen_authorities"]["config"]["path"], label="native config")
    constructor = config["constructor"]
    torch.manual_seed(12031986)
    model = scooby.Scooby(
        BorzoiConfig(**native_config),
        cell_emb_dim=constructor["cell_emb_dim"],
        embedding_dim=constructor["embedding_dim"],
        n_tracks=constructor["n_tracks"],
        disable_cache=constructor["disable_cache"],
        use_transform_borzoi_emb=constructor["use_transform_borzoi_emb"],
    )
    model.eval()
    state = model.state_dict()
    header = load_json(root / config["frozen_authorities"]["header"]["path"], label="bounded header")
    expected = {
        row["name"]: (tuple(row["shape"]), row["dtype"], row["parameter_count"])
        for row in header["tensors"]
    }
    dtype_names = {torch.float32: "F32", torch.int64: "I64"}
    observed = {
        name: (tuple(tensor.shape), dtype_names.get(tensor.dtype, str(tensor.dtype)), tensor.numel())
        for name, tensor in state.items()
    }
    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    mismatched = sorted(name for name in set(expected) & set(observed) if expected[name] != observed[name])
    require(not missing and not unexpected and not mismatched, "runtime state schema differs from checkpoint header")
    require(len(observed) == constructor["expected_state_dict_tensors"], "state tensor count differs")
    require(sum(value[2] for value in observed.values()) == constructor["expected_state_dict_elements_including_buffers"], "state element count differs")

    cell_context = torch.linspace(-1.0, 1.0, steps=constructor["cell_emb_dim"], dtype=torch.float32).reshape(1, 1, -1)
    with torch.inference_mode():
        weights_a, biases_a = model.forward_cell_embs_only(cell_context)
        weights_b, biases_b = model.forward_cell_embs_only(cell_context)
        sequence_embedding = torch.linspace(-0.25, 0.25, steps=constructor["embedding_dim"] * 8, dtype=torch.float32).reshape(1, constructor["embedding_dim"], 8)
        profile_a = model.forward_convs_on_emb(sequence_embedding, weights_a, biases_a)
        profile_b = model.forward_convs_on_emb(sequence_embedding, weights_b, biases_b)
    require(torch.equal(weights_a, weights_b) and torch.equal(biases_a, biases_b), "decoder repeat differs")
    require(torch.equal(profile_a, profile_b), "profile repeat differs")
    require(tuple(weights_a.shape) == (1, 3, 1920, 1), "decoder weight shape differs")
    require(tuple(biases_a.shape) == (1, 3), "decoder bias shape differs")
    require(tuple(profile_a.shape) == (1, 8, 3), "synthetic profile shape differs")
    require(torch.isfinite(profile_a).all().item() and torch.all(profile_a > 0).item(), "softplus profile contract differs")
    prediction = profile_a.detach().cpu().numpy().astype("<f4", copy=False)
    return {
        "runtime_versions": versions,
        "state_dict_tensors": len(observed),
        "state_dict_elements_including_buffers": sum(value[2] for value in observed.values()),
        "schema_missing": missing,
        "schema_unexpected": unexpected,
        "schema_mismatched": mismatched,
        "random_weight_decoder_probe_shape": list(prediction.shape),
        "random_weight_decoder_probe_sha256": hashlib.sha256(prediction.tobytes(order="C")).hexdigest(),
        "deterministic_repeat_bit_identical": True,
        "checkpoint_weights_loaded": False,
        "full_sequence_forward_executed": False,
        "peft_training_or_mutation_called": False,
    }


def audit(root: Path, config_path: Path, source_dir: Path, output: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    config = load_json(config_path.resolve(strict=True), label="schema gate config")
    require(config.get("schema_version") == SCHEMA, "schema gate identity differs")
    require(not output.exists() and not output.is_symlink(), "refusing to overwrite output")
    policy = config["gate_policy"]
    require(
        not policy["project_data_read_allowed"]
        and not policy["sealed_data_read_allowed"]
        and not policy["checkpoint_payload_downloaded_in_this_gate"]
        and not policy["checkpoint_deserialized_in_this_gate"]
        and not policy["full_sequence_forward_allowed"]
        and policy["random_weight_decoder_probe_allowed"],
        "schema gate firewall differs",
    )
    authorities = validate_authorities(root, config)
    sources = validate_sources(source_dir.resolve(strict=True), config)
    runtime = runtime_probe(root, source_dir, config)
    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": "masld-bench-scooby-epicardioids-runtime-schema-receipt-v1",
        "status": "pass_checkpoint_download_and_strict_restore_probe_allowed",
        "config_sha256": digest(config_path),
        "frozen_authorities": authorities,
        "official_sources": sources,
        "runtime_probe": runtime,
        "task_scope": config["task_scope"],
        "checkpoint_payload_downloaded": False,
        "checkpoint_deserialized": False,
        "checkpoint_forward_executed": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "rna_conditioned_atac_eligible": False,
        "observed_multiome_development_comparator_in_principle": True,
        "open_champion_eligible": False,
        "next_gate": policy["next_if_pass"],
    }
    write_json_exclusive(output / "receipt.json", receipt, mode=0o640)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipt = audit(args.root, args.config, args.source_dir, args.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
