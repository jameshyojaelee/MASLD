#!/usr/bin/env python3
"""Inspect an exact PyTorch checkpoint through the restricted weights-only loader."""

from __future__ import annotations

import argparse
from collections import Counter
from collections.abc import Mapping
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


class RestrictedInspectionError(ValueError):
    """Raised when a checkpoint or loaded value does not meet the inspection requirements."""


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _type_name(value: Any) -> str:
    value_type = type(value)
    return f"{value_type.__module__}.{value_type.__qualname__}"


def _scalar_record(value: Any) -> dict[str, Any]:
    if value is None or isinstance(value, (bool, int, float)):
        return {"type": _type_name(value), "value": value}
    if isinstance(value, str):
        encoded = value.encode("utf-8")
        return {
            "type": "builtins.str",
            "length": len(value),
            "sha256": sha256(encoded).hexdigest(),
            "value": value if len(value) <= 160 else None,
        }
    if isinstance(value, bytes):
        return {
            "type": "builtins.bytes",
            "length": len(value),
            "sha256": sha256(value).hexdigest(),
        }
    raise RestrictedInspectionError(f"unsupported scalar type: {_type_name(value)}")


def inspect_tree(payload: Any, *, maximum_nodes: int = 2_000_000) -> dict[str, Any]:
    import torch

    tensors: list[dict[str, Any]] = []
    containers: list[dict[str, Any]] = []
    scalars: list[dict[str, Any]] = []
    type_counts: Counter[str] = Counter()
    nodes_seen = 0

    def visit(value: Any, path: str, depth: int) -> None:
        nonlocal nodes_seen
        nodes_seen += 1
        if nodes_seen > maximum_nodes or depth > 64:
            raise RestrictedInspectionError("checkpoint object graph exceeds inspection limits")
        type_counts[_type_name(value)] += 1

        if isinstance(value, torch.Tensor):
            tensors.append(
                {
                    "path": path,
                    "shape": list(value.shape),
                    "dtype": str(value.dtype),
                    "device": str(value.device),
                    "requires_grad": bool(value.requires_grad),
                    "numel": int(value.numel()),
                    "storage_bytes": int(value.numel() * value.element_size()),
                }
            )
            return
        if isinstance(value, Mapping):
            if any(not isinstance(key, (str, int)) for key in value):
                raise RestrictedInspectionError(
                    f"mapping at {path} contains a non-string/non-integer key"
                )
            ordered_keys = sorted(value, key=lambda key: (str(type(key)), str(key)))
            containers.append(
                {
                    "path": path,
                    "type": _type_name(value),
                    "length": len(value),
                    "keys": [str(key) for key in ordered_keys],
                }
            )
            for key in ordered_keys:
                visit(value[key], f"{path}.{key}", depth + 1)
            return
        if isinstance(value, (list, tuple)):
            containers.append(
                {
                    "path": path,
                    "type": _type_name(value),
                    "length": len(value),
                    "item_type_counts": dict(
                        sorted(Counter(_type_name(item) for item in value).items())
                    ),
                }
            )
            if len(value) <= 256:
                for index, item in enumerate(value):
                    visit(item, f"{path}[{index}]", depth + 1)
            return
        if isinstance(value, (type(None), bool, int, float, str, bytes)):
            record = _scalar_record(value)
            record["path"] = path
            scalars.append(record)
            return
        if isinstance(value, (torch.Size, torch.dtype, torch.device)):
            scalars.append(
                {"path": path, "type": _type_name(value), "value": str(value)}
            )
            return
        raise RestrictedInspectionError(
            f"weights-only loader returned unsupported type at {path}: {_type_name(value)}"
        )

    visit(payload, "checkpoint", 0)
    return {
        "node_count": nodes_seen,
        "type_counts": dict(sorted(type_counts.items())),
        "tensor_count": len(tensors),
        "tensor_numel": sum(row["numel"] for row in tensors),
        "tensor_storage_bytes": sum(row["storage_bytes"] for row in tensors),
        "tensors": tensors,
        "containers": containers,
        "scalars": scalars,
    }


def inspect(
    checkpoint: Path,
    output: Path,
    *,
    expected_sha256: str,
    expected_size: int,
) -> dict[str, Any]:
    if output.exists() or checkpoint.is_symlink() or not checkpoint.is_file():
        raise RestrictedInspectionError("checkpoint request differs")
    if checkpoint.stat().st_size != expected_size:
        raise RestrictedInspectionError("checkpoint size differs")
    observed_sha256 = _sha256_file(checkpoint)
    if observed_sha256 != expected_sha256:
        raise RestrictedInspectionError("checkpoint SHA-256 differs")

    import torch

    payload = torch.load(
        checkpoint,
        map_location="cpu",
        weights_only=True,
        mmap=True,
    )
    tree = inspect_tree(payload)
    result = {
        "schema_version": "masld-bench-torch-weights-only-inspection-v1",
        "status": "pass",
        "checkpoint": checkpoint.name,
        "sha256": observed_sha256,
        "size_bytes": expected_size,
        "torch_version": torch.__version__,
        "weights_only": True,
        "map_location": "cpu",
        "mmap": True,
        **tree,
    }
    output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--expected-size", type=int, required=True)
    arguments = parser.parse_args()
    result = inspect(
        arguments.checkpoint,
        arguments.output,
        expected_sha256=arguments.expected_sha256,
        expected_size=arguments.expected_size,
    )
    print(
        json.dumps(
            {
                key: value
                for key, value in result.items()
                if key not in {"tensors", "containers", "scalars", "type_counts"}
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
