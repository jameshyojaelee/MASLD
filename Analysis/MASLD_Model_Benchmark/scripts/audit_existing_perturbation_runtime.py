#!/usr/bin/env python3
"""Audit one existing perturbation-model environment without installing packages."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
from typing import Any


LANES = {
    "gears": {"modules": ["gears", "torch", "torch_geometric", "anndata"], "distributions": ["cell-gears", "torch", "torch-geometric", "anndata"]},
    "scgpt": {"modules": ["scgpt", "torch", "anndata"], "distributions": ["scgpt", "torch", "anndata"]},
    "scgen": {"modules": ["scgen", "torch", "scvi", "anndata"], "distributions": ["scgen", "torch", "scvi-tools", "anndata"]},
    "cpa": {"modules": ["cpa", "torch", "scvi", "anndata"], "distributions": ["cpa-tools", "torch", "scvi-tools", "anndata"]},
    "cellot": {"modules": ["cellot", "jax", "ott", "anndata"], "distributions": ["cellot", "jax", "ott-jax", "anndata"]},
    "genepert": {"modules": ["GenePert", "torch", "anndata"], "distributions": ["GenePert", "torch", "anndata"]},
    "perturblib": {"modules": ["perturb_lib", "torch", "anndata"], "distributions": ["perturblib", "torch", "anndata"]},
}


def file_digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def distribution(name: str) -> dict[str, Any]:
    try:
        item = importlib.metadata.distribution(name)
    except importlib.metadata.PackageNotFoundError:
        return {"name": name, "installed": False}
    root = Path(item._path)
    record = root / "RECORD"
    return {
        "name": name,
        "installed": True,
        "version": item.version,
        "metadata_path": root.as_posix(),
        "record_sha256": file_digest(record) if record.is_file() else None,
    }


def module(name: str) -> dict[str, Any]:
    try:
        item = importlib.import_module(name)
    except Exception as error:
        return {"name": name, "imported": False, "error_type": type(error).__name__, "error": str(error)[:2000]}
    file_value = getattr(item, "__file__", None)
    path = Path(file_value).resolve() if file_value else None
    result: dict[str, Any] = {
        "name": name,
        "imported": True,
        "version": str(getattr(item, "__version__", "UNDECLARED")),
        "path": path.as_posix() if path else None,
        "file_sha256": file_digest(path) if path and path.is_file() else None,
    }
    if name == "torch":
        result.update(
            {
                "cuda_available": bool(item.cuda.is_available()),
                "cuda_runtime": str(item.version.cuda),
                "device_count": int(item.cuda.device_count()),
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lane", choices=tuple(LANES), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError("runtime audit would overwrite output")
    spec = LANES[args.lane]
    distributions = [distribution(name) for name in spec["distributions"]]
    modules = [module(name) for name in spec["modules"]]
    primary_imported = bool(modules and modules[0].get("imported"))
    result = {
        "schema_version": "masld-bench-existing-perturbation-runtime-audit-v1",
        "lane": args.lane,
        "status": "observed_not_admitted",
        "primary_imported": primary_imported,
        "python": {"executable": str(Path(sys.executable).resolve()), "version": sys.version, "platform": platform.platform()},
        "distributions": distributions,
        "modules": modules,
        "installation_performed": False,
        "forward_fixture_performed": False,
        "admission_note": "A mutable environment import census does not admit a faithful runtime or model forward pass.",
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "runtime_audit.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"lane": args.lane, "primary_imported": primary_imported}, sort_keys=True))


if __name__ == "__main__":
    main()
