#!/usr/bin/env python3
"""Acquire only the hash-pinned LS-GKM build surface, never historical data."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import ssl
import tomllib
from typing import Any, Mapping
from urllib.parse import urlparse
from urllib.request import Request, urlopen


SCHEMA = "masld-bench-lsgkm-runtime-fixture-v1"
MAX_SOURCE_FILE_BYTES = 4 * 1024 * 1024


class SourceAcquisitionError(RuntimeError):
    """Raised when remote identity or the included source surface differs."""


def load_config(path: Path) -> dict[str, Any]:
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as error:
        raise SourceAcquisitionError(f"invalid runtime config: {error}") from error
    if not isinstance(value, dict) or value.get("schema_version") != SCHEMA:
        raise SourceAcquisitionError("runtime config schema differs")
    return value


def canonical_sha256(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def admitted_records(config: Mapping[str, Any]) -> list[dict[str, str]]:
    source = config.get("source", {})
    admitted = source.get("admitted_paths", [])
    records = source.get("files", [])
    if (
        not isinstance(admitted, list)
        or not isinstance(records, list)
        or len(admitted) != len(records)
        or set(admitted) != {record.get("path") for record in records}
    ):
        raise SourceAcquisitionError("source path census differs")
    prohibited = tuple(source.get("prohibited_path_prefixes", ()))
    normalized: list[dict[str, str]] = []
    for record in records:
        relative = record.get("path")
        expected = record.get("sha256")
        path = PurePosixPath(str(relative))
        if (
            path.is_absolute()
            or ".." in path.parts
            or str(path) != relative
            or any(relative.startswith(prefix) for prefix in prohibited)
            or not isinstance(expected, str)
            or len(expected) != 64
        ):
            raise SourceAcquisitionError("unsafe or prohibited source record")
        normalized.append({"path": relative, "sha256": expected})
    return normalized


def fetch_bytes(url: str) -> bytes:
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != "raw.githubusercontent.com":
        raise SourceAcquisitionError("source URL host or scheme differs")
    request = Request(url, headers={"User-Agent": "masld-bench-pinned-source/1"})
    context = ssl.create_default_context()
    try:
        with urlopen(request, timeout=120, context=context) as response:
            if response.geturl() != url or response.status != 200:
                raise SourceAcquisitionError("source redirect or status differs")
            length = response.headers.get("Content-Length")
            if length is not None and int(length) > MAX_SOURCE_FILE_BYTES:
                raise SourceAcquisitionError("source file exceeds byte bound")
            payload = response.read(MAX_SOURCE_FILE_BYTES + 1)
    except SourceAcquisitionError:
        raise
    except Exception as error:
        raise SourceAcquisitionError(f"source acquisition failed: {error}") from error
    if not payload or len(payload) > MAX_SOURCE_FILE_BYTES:
        raise SourceAcquisitionError("source file is empty or exceeds byte bound")
    return payload


def acquire(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise SourceAcquisitionError("source output exists")
    config = load_config(arguments.config)
    if (
        config.get("historical_weights_allowed") is not False
        or config.get("historical_implementation_allowed") is not False
        or config.get("outcome_access_authorized") is not False
        or config.get("biological_data_access_authorized") is not False
    ):
        raise SourceAcquisitionError("source or outcome authorization differs")
    base_url = config["source"]["raw_base_url"].rstrip("/")
    records = admitted_records(config)
    arguments.output.mkdir(parents=True, mode=0o750)
    inventory: list[dict[str, Any]] = []
    for record in records:
        relative = record["path"]
        url = f"{base_url}/{relative}"
        payload = fetch_bytes(url)
        observed = canonical_sha256(payload)
        if observed != record["sha256"]:
            raise SourceAcquisitionError(f"source hash differs: {relative}")
        target = arguments.output / relative
        target.parent.mkdir(parents=True, exist_ok=True, mode=0o750)
        with target.open("xb") as handle:
            handle.write(payload)
        inventory.append(
            {
                "path": relative,
                "url": url,
                "sha256": observed,
                "size_bytes": len(payload),
            }
        )
    receipt = {
        "schema_version": "masld-bench-lsgkm-pinned-source-receipt-v1",
        "status": "pass_exact_source_only_no_historical_weights",
        "repository": "Dongwon-Lee/lsgkm",
        "revision": config["source_revision"],
        "revision_tree_sha256": config["source_tree_sha256"],
        "files": inventory,
        "admitted_file_count": len(inventory),
        "historical_weights_acquired": False,
        "historical_implementation_acquired": False,
        "repository_cloned": False,
        "tests_directory_acquired": False,
        "outcomes_read": False,
        "biological_data_read": False,
    }
    with (arguments.output / "source_receipt.json").open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(acquire(parse_args()), sort_keys=True))
