#!/usr/bin/env python3
"""Audit an existing cell-foundation environment without installing packages."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
from typing import Any


LANES = {
    "decima": {
        "distributions": ("torch", "transformers", "anndata", "scanpy", "numpy", "scikit-learn"),
        "modules": ("torch", "transformers", "anndata", "scanpy", "numpy", "sklearn"),
        "intended_models": ("geneformer", "scgpt", "uce"),
    },
    "scimilarity": {
        "distributions": ("torch", "scimilarity", "anndata", "scanpy", "numpy"),
        "modules": ("torch", "scimilarity", "anndata", "scanpy", "numpy"),
        "intended_models": ("scimilarity",),
    },
    "scprint310": {
        "distributions": ("torch", "scprint", "anndata", "scanpy", "numpy"),
        "modules": ("torch", "scprint", "anndata", "scanpy", "numpy"),
        "intended_models": ("scprint_v1_5_medium",),
    },
    "scprint2x": {
        "distributions": ("torch", "scprint", "scprint2", "anndata", "scanpy", "numpy"),
        "modules": ("torch", "scprint", "scprint2", "anndata", "scanpy", "numpy"),
        "intended_models": ("scprint2_small_v2",),
    },
}


class RuntimeAuditError(ValueError):
    """Raised when the read-only audit request differs."""


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _distribution(name: str) -> dict[str, Any]:
    try:
        distribution = importlib.metadata.distribution(name)
    except importlib.metadata.PackageNotFoundError:
        return {"name": name, "installed": False}
    metadata_path = Path(distribution._path)  # read-only dist-info identity
    record = metadata_path / "RECORD"
    return {
        "name": name,
        "installed": True,
        "version": distribution.version,
        "metadata_path": metadata_path.as_posix(),
        "record_sha256": _sha256_file(record) if record.is_file() else None,
    }


def _module(name: str) -> dict[str, Any]:
    try:
        module = importlib.import_module(name)
    except Exception as error:  # preserve dependency/import failures as evidence
        return {
            "name": name,
            "imported": False,
            "error_type": type(error).__name__,
            "error": str(error)[:1000],
        }
    path_value = getattr(module, "__file__", None)
    path = Path(path_value).resolve() if path_value else None
    result: dict[str, Any] = {
        "name": name,
        "imported": True,
        "module_version": str(getattr(module, "__version__", "UNDECLARED")),
        "module_path": path.as_posix() if path else None,
        "module_file_sha256": _sha256_file(path) if path and path.is_file() else None,
    }
    if name == "torch":
        cuda_available = bool(module.cuda.is_available())
        result["cuda_available"] = cuda_available
        result["cuda_runtime"] = str(module.version.cuda)
        result["cudnn_version"] = module.backends.cudnn.version()
        result["device_count_reported"] = int(module.cuda.device_count())
        result["devices"] = []
        if cuda_available:
            result["devices"] = [
                {
                    "index": index,
                    "name": module.cuda.get_device_name(index),
                    "capability": list(module.cuda.get_device_capability(index)),
                }
                for index in range(module.cuda.device_count())
            ]
    return result


def audit(lane: str, output: Path) -> dict[str, Any]:
    if lane not in LANES or output.exists():
        raise RuntimeAuditError("runtime audit request differs")
    spec = LANES[lane]
    distributions = [_distribution(name) for name in spec["distributions"]]
    modules = [_module(name) for name in spec["modules"]]
    torch_result = next(item for item in modules if item["name"] == "torch")
    result = {
        "schema_version": "masld-bench-existing-cell-foundation-runtime-audit-v1",
        "status": "observed_not_admitted",
        "lane": lane,
        "intended_models": list(spec["intended_models"]),
        "network_policy": {
            "HF_HUB_OFFLINE": os.environ.get("HF_HUB_OFFLINE"),
            "TRANSFORMERS_OFFLINE": os.environ.get("TRANSFORMERS_OFFLINE"),
            "WANDB_MODE": os.environ.get("WANDB_MODE"),
        },
        "python": {
            "version": sys.version,
            "executable": Path(sys.executable).resolve().as_posix(),
            "platform": platform.platform(),
        },
        "distributions": distributions,
        "modules": modules,
        "gpu_runtime_observed": bool(
            torch_result.get("imported") and torch_result.get("cuda_available")
        ),
        "admission_note": (
            "This is a read-only census of a mutable pre-existing environment. "
            "It is not a faithful-runtime or model-forward parity admission."
        ),
    }
    output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lane", choices=tuple(LANES), required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = audit(arguments.lane, arguments.output)
    print(json.dumps({key: result[key] for key in ("status", "lane", "gpu_runtime_observed")}))


if __name__ == "__main__":
    main()
