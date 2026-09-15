#!/usr/bin/env python3
"""Convert included PyTorch ZIP checkpoints using the weights-only unpickler."""

from __future__ import annotations

import argparse
from collections.abc import Mapping
from hashlib import sha256
import io
import json
from pathlib import Path
from typing import Any

import torch
from safetensors.torch import save_file


class WeightsOnlyConversionError(ValueError):
    """Raised when a checkpoint cannot be reduced to one tensor state dict."""


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _object_schema(value: Any, depth: int = 0) -> object:
    if torch.is_tensor(value):
        return {
            "type": "tensor",
            "shape": list(value.shape),
            "dtype": str(value.dtype),
        }
    if isinstance(value, io.BytesIO):
        payload = value.getbuffer()
        return {
            "type": "BytesIO",
            "size_bytes": len(payload),
            "sha256": sha256(payload).hexdigest(),
        }
    if depth >= 3:
        return {"type": type(value).__name__, "depth_truncated": True}
    if isinstance(value, Mapping):
        return {
            "type": type(value).__name__,
            "length": len(value),
            "items": {
                str(key): _object_schema(item, depth + 1)
                for key, item in list(value.items())[:20]
            },
            "items_truncated": len(value) > 20,
        }
    if isinstance(value, (list, tuple)):
        return {
            "type": type(value).__name__,
            "length": len(value),
            "items": [_object_schema(item, depth + 1) for item in value[:20]],
            "items_truncated": len(value) > 20,
        }
    if value is None or isinstance(value, (str, int, float, bool)):
        return {"type": type(value).__name__}
    raise WeightsOnlyConversionError(
        f"weights-only checkpoint contains an unsupported object: {type(value)!r}"
    )


def _tensor_mapping(value: Any) -> dict[str, torch.Tensor] | None:
    if not isinstance(value, Mapping) or not value:
        return None
    if not all(isinstance(key, str) and torch.is_tensor(item) for key, item in value.items()):
        return None
    return dict(value)


def _tensor_mapping_with_safe_metadata(
    value: Any,
) -> tuple[dict[str, torch.Tensor], list[dict[str, object]]] | None:
    if not isinstance(value, Mapping) or not value:
        return None
    tensors: dict[str, torch.Tensor] = {}
    metadata: list[dict[str, object]] = []
    for key, item in value.items():
        if not isinstance(key, str):
            return None
        if torch.is_tensor(item):
            tensors[key] = item
        elif isinstance(item, io.BytesIO):
            payload = item.getbuffer()
            metadata.append(
                {
                    "key": key,
                    "type": "BytesIO",
                    "size_bytes": len(payload),
                    "sha256": sha256(payload).hexdigest(),
                }
            )
        elif item is None or isinstance(item, (str, int, float, bool)):
            metadata.append({"key": key, "type": type(item).__name__})
        else:
            return None
    return (tensors, metadata) if tensors else None


def _select_state_dict(
    value: Any,
) -> tuple[dict[str, torch.Tensor], str, list[dict[str, object]]]:
    direct = _tensor_mapping(value)
    if direct is not None:
        return direct, "root", []
    if not isinstance(value, Mapping):
        raise WeightsOnlyConversionError("checkpoint root is not a mapping")
    candidates: list[
        tuple[str, dict[str, torch.Tensor], list[dict[str, object]]]
    ] = []

    def visit(item: Any, path: str, depth: int) -> None:
        selected = _tensor_mapping_with_safe_metadata(item)
        if selected is not None:
            tensors, metadata = selected
            candidates.append((path, tensors, metadata))
        if depth >= 5 or not isinstance(item, Mapping):
            return
        for key, nested in item.items():
            if isinstance(key, str) and isinstance(nested, Mapping):
                visit(nested, f"{path}.{key}" if path else key, depth + 1)

    visit(value, "root", 0)
    candidates.sort(
        key=lambda candidate: sum(tensor.numel() for tensor in candidate[1].values()),
        reverse=True,
    )
    if not candidates:
        raise WeightsOnlyConversionError(
            "checkpoint has no tensor mapping with only inert metadata"
        )
    if len(candidates) > 1:
        largest = sum(tensor.numel() for tensor in candidates[0][1].values())
        second = sum(tensor.numel() for tensor in candidates[1][1].values())
        if largest <= 2 * second:
            census = [
                {
                    "path": path,
                    "tensor_keys": len(tensors),
                    "tensor_elements": sum(tensor.numel() for tensor in tensors.values()),
                }
                for path, tensors, _ in candidates[:20]
            ]
            raise WeightsOnlyConversionError(
                f"checkpoint tensor mapping is ambiguous: {census!r}"
            )
    path, tensors, metadata = candidates[0]
    return tensors, path, metadata


def _storage_identity(tensor: torch.Tensor) -> tuple[int, int]:
    storage = tensor.untyped_storage()
    return storage.data_ptr(), storage.nbytes()


def _deduplicate_aliases(
    state: dict[str, torch.Tensor],
) -> tuple[dict[str, torch.Tensor], list[dict[str, object]]]:
    groups: dict[tuple[int, int], list[str]] = {}
    for key, tensor in state.items():
        groups.setdefault(_storage_identity(tensor), []).append(key)
    output: dict[str, torch.Tensor] = {}
    aliases: list[dict[str, object]] = []
    for storage_identity, keys in groups.items():
        views: dict[tuple[object, ...], list[str]] = {}
        for key in sorted(keys):
            tensor = state[key]
            signature: tuple[object, ...] = (
                tuple(tensor.shape),
                tuple(tensor.stride()),
                tensor.storage_offset(),
                str(tensor.dtype),
            )
            views.setdefault(signature, []).append(key)
        materialized_views: list[dict[str, object]] = []
        for signature, view_keys in views.items():
            canonical_key = view_keys[0]
            # Safetensors forbids shared storage. Preserve every distinct view
            # as an independent contiguous tensor while omitting only names
            # that point to the exact same source view.
            output[canonical_key] = state[canonical_key].contiguous().clone()
            materialized_views.append(
                {
                    "canonical": canonical_key,
                    "omitted_identical_view_aliases": view_keys[1:],
                    "shape": list(signature[0]),
                    "stride": list(signature[1]),
                    "storage_offset": signature[2],
                    "dtype": signature[3],
                }
            )
        if len(keys) > 1:
            aliases.append(
                {
                    "source_storage_nbytes": storage_identity[1],
                    "source_keys": sorted(keys),
                    "nonidentical_views_materialized_separately": len(views) > 1,
                    "materialized_views": materialized_views,
                }
            )
    return output, aliases


def convert(
    source: Path,
    output: Path,
    expected_sha256: str,
    expected_size: int,
    model_id: str,
    source_container: str = "PyTorch_ZIP_with_pickle_metadata",
    use_mmap: bool = True,
) -> dict[str, object]:
    if output.exists() or source.is_symlink() or not source.is_file():
        raise WeightsOnlyConversionError("weights-only conversion request differs")
    if source.stat().st_size != expected_size or _sha256_file(source) != expected_sha256:
        raise WeightsOnlyConversionError("source checkpoint identity differs")
    output.mkdir(mode=0o750)
    safe_globals: list[str] = []
    if model_id == "evo":
        # Evo 2 wraps its tensor mapping in an in-memory byte stream. BytesIO
        # carries bytes and has no external side effects; no other class is
        # included if the weights-only unpickler encounters it.
        torch.serialization.add_safe_globals([io.BytesIO])
        safe_globals.append("_io.BytesIO")
    loaded = torch.load(
        source,
        map_location="cpu",
        weights_only=True,
        mmap=use_mmap,
    )
    schema = _object_schema(loaded)
    (output / "source_object_schema.json").write_text(
        json.dumps(schema, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    state, selected_path, ignored_metadata = _select_state_dict(loaded)
    if not state or any(not key or "\x00" in key for key in state):
        raise WeightsOnlyConversionError("state-dict key census differs")
    converted, aliases = _deduplicate_aliases(state)
    if not all(tensor.device.type == "cpu" for tensor in converted.values()):
        raise WeightsOnlyConversionError("state dict is not CPU resident")
    dtype_counts: dict[str, int] = {}
    element_count = 0
    for tensor in converted.values():
        dtype_counts[str(tensor.dtype)] = dtype_counts.get(str(tensor.dtype), 0) + 1
        element_count += tensor.numel()
    target = output / f"{model_id}.model.safetensors"
    save_file(
        converted,
        str(target),
        metadata={
            "format": "pt",
            "source_checkpoint_sha256": expected_sha256,
            "source_state_dict_path": selected_path,
        },
    )
    receipt = {
        "schema_version": "masld-bench-dna-lm-weights-only-conversion-v1",
        "status": "pass",
        "model_id": model_id,
        "source_checkpoint_sha256": expected_sha256,
        "source_checkpoint_size_bytes": expected_size,
        "source_container": source_container,
        "loader": f"torch.load_weights_only_true_mmap_{str(use_mmap).lower()}",
        "explicit_weights_only_safe_globals": safe_globals,
        "source_object_schema": schema,
        "selected_state_dict_path": selected_path,
        "ignored_inert_state_dict_metadata": ignored_metadata,
        "source_tensor_key_count": len(state),
        "converted_tensor_key_count": len(converted),
        "converted_tensor_element_count": element_count,
        "converted_tensor_dtype_counts": dtype_counts,
        "shared_storage_aliases": aliases,
        "converted_safetensors_sha256": _sha256_file(target),
        "converted_safetensors_size_bytes": target.stat().st_size,
        "custom_model_code_imported": False,
        "model_instantiated": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "runtime_network_allowed": False,
    }
    (output / "conversion_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--expected-size", type=int, required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument(
        "--source-container",
        choices=("PyTorch_ZIP_with_pickle_metadata", "PyTorch_legacy_pickle_stream"),
        default="PyTorch_ZIP_with_pickle_metadata",
    )
    parser.add_argument("--disable-mmap", action="store_true")
    args = parser.parse_args()
    convert(
        args.source,
        args.output,
        args.expected_sha256,
        args.expected_size,
        args.model_id,
        args.source_container,
        not args.disable_mmap,
    )


if __name__ == "__main__":
    main()
