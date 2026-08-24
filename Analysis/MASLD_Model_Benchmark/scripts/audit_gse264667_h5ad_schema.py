#!/usr/bin/env python3
"""Audit GSE264667 HepG2 H5AD topology without reading dataset values."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import h5py


class H5adAuditError(RuntimeError):
    """Raised when the H5AD topology violates bounded safe-loading rules."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def attribute_schema(obj: h5py.Group | h5py.Dataset) -> list[dict[str, Any]]:
    records = []
    for name in sorted(obj.attrs):
        attribute = obj.attrs.get_id(name)
        records.append(
            {
                "name": name,
                "shape": list(attribute.shape),
                "dtype": str(attribute.dtype),
            }
        )
    return records


def audit(source: Path) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    visited_groups: set[int] = set()

    def walk(group: h5py.Group, prefix: str) -> None:
        address = int(h5py.h5o.get_info(group.id).addr)
        if address in visited_groups:
            raise H5adAuditError(f"hard-linked group cycle or alias at {prefix}")
        visited_groups.add(address)
        for key in sorted(group.keys()):
            path = f"{prefix}/{key}" if prefix else f"/{key}"
            link = group.get(key, getlink=True)
            if not isinstance(link, h5py.HardLink):
                raise H5adAuditError(f"non-hard HDF5 link prohibited: {path}")
            obj = group.get(key)
            if isinstance(obj, h5py.Group):
                records.append(
                    {
                        "path": path,
                        "kind": "group",
                        "attributes": attribute_schema(obj),
                    }
                )
                walk(obj, path)
            elif isinstance(obj, h5py.Dataset):
                records.append(
                    {
                        "path": path,
                        "kind": "dataset",
                        "shape": list(obj.shape),
                        "dtype": str(obj.dtype),
                        "chunks": list(obj.chunks) if obj.chunks else None,
                        "compression": obj.compression,
                        "stored_size_bytes": int(obj.id.get_storage_size()),
                        "attributes": attribute_schema(obj),
                    }
                )
            else:
                raise H5adAuditError(f"unsupported HDF5 object: {path}")
            if len(records) > 50_000:
                raise H5adAuditError("HDF5 object count exceeds 50,000 ceiling")

    with h5py.File(source, mode="r", locking=True) as handle:
        root_attributes = attribute_schema(handle)
        walk(handle, "")
    matrix_paths = [
        record for record in records
        if record["path"] in {"/X", "/X/data", "/raw/X", "/raw/X/data"}
    ]
    return {
        "schema_version": "masld-bench-gse264667-h5ad-schema-audit-v1",
        "source_path": source.as_posix(),
        "source_sha256": digest(source),
        "source_size_bytes": source.stat().st_size,
        "h5py_version": h5py.__version__,
        "root_attributes": root_attributes,
        "object_count": len(records),
        "objects": records,
        "matrix_storage_records": matrix_paths,
        "dataset_values_read": False,
        "expression_values_read": False,
        "external_or_soft_links_allowed": False,
        "biological_unit_status": "single_HepG2_transduced_pool_GEM_groups_not_replicates",
        "gse313774_accessed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.source.is_file() or args.output.exists():
        raise H5adAuditError("H5AD schema audit request differs")
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "h5ad_schema.json").write_text(
        json.dumps(audit(args.source), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
