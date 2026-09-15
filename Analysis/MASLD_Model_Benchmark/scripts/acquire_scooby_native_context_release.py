#!/usr/bin/env python3
"""Selectively acquire released Scooby native-context rows from a remote ZIP."""

from __future__ import annotations

import argparse
import binascii
import hashlib
import json
from pathlib import Path, PurePath
import struct
from typing import Any, Mapping
import zlib

from masld_bench.artifacts import ArtifactError, verify_frozen_tree, write_bytes_exclusive, write_json_exclusive
from scripts.scooby_remote_zip_inventory import HTTPRangeReader


SCHEMA = "masld-bench-scooby-epicardioids-native-context-release-v1"
LOCAL_SIGNATURE = b"PK\x03\x04"


class ScoobyNativeContextReleaseError(RuntimeError):
    """Raised when a released context member or its boundary differs."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ScoobyNativeContextReleaseError(message)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    require(path.is_file() and not path.is_symlink(), f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScoobyNativeContextReleaseError(f"invalid {label}: {error}") from error
    require(isinstance(value, dict), f"{label} is not an object")
    return value


def validate_boundary(config: Mapping[str, Any]) -> None:
    boundary = config["release_boundary"]
    require(boundary["native_context_rows_released"], "released native context rows were not admitted")
    require(not boundary["exact_scglue_encoder_weights_released_in_archive"], "unreleased scGLUE encoder was asserted")
    require(not boundary["exact_guidance_graph_released_in_archive"], "unreleased guidance graph was asserted")
    require(not boundary["frozen_query_transform_released_in_archive"], "unreleased query transform was asserted")
    require(not boundary["coordinate_compatible_liver_query_supported"], "coordinate-compatible liver query was asserted")
    require(not boundary["rna_conditioned_atac_eligible"], "observed-ATAC context cannot enter RNA-conditioned lane")
    require(boundary["observed_multiome_context_only"], "observed-multiome task boundary differs")
    require(not boundary["project_data_read_allowed"] and not boundary["sealed_data_read_allowed"], "data firewall differs")
    require(not boundary["full_archive_download_allowed"], "selective acquisition cannot download the full archive")


def validate_inventory(config: Mapping[str, Any], inventory: Mapping[str, Any]) -> list[dict[str, Any]]:
    archive = config["archive"]
    require(inventory["remote_url"] == archive["url"], "archive URL differs")
    require(inventory["archive_size_bytes"] == archive["size_bytes"], "archive size differs")
    require(inventory["declared_archive_md5"] == archive["declared_md5"], "archive MD5 declaration differs")
    require(not inventory["declared_archive_md5_verified_from_full_payload"], "unexpected full-archive claim")
    require(inventory["central_directory_offset"] == archive["central_directory_offset"], "central directory offset differs")
    require(inventory["central_directory_size_bytes"] == archive["central_directory_size_bytes"], "central directory size differs")
    require(inventory["entry_count"] == archive["entry_count"], "archive entry count differs")
    indexed = {member["path"]: member for member in inventory["members"]}
    for expected in config["members"]:
        observed = indexed.get(expected["path"])
        require(observed is not None, f"released member is absent: {expected['path']}")
        for key in ("local_header_offset", "compression", "compressed_size", "uncompressed_size", "crc32"):
            require(observed[key] == expected[key], f"released member contract differs: {expected['path']} {key}")
    prefix = "training_data/epicardioids_training_data/"
    epicardioid_members = [member for member in inventory["members"] if member["path"].startswith(prefix) and not member["directory"]]
    encoder_candidates = [
        member for member in epicardioid_members
        if any(token in member["path"].lower() for token in ("scglue", "guidance", "query_model", "reference_model"))
        or PurePath(member["path"]).suffix.lower() in {".pt", ".ckpt", ".pth"}
    ]
    require(not encoder_candidates, "an encoder candidate exists and needs separate adjudication")
    return epicardioid_members


def extract_member(reader: Any, member: Mapping[str, Any]) -> tuple[bytes, dict[str, Any]]:
    offset = int(member["local_header_offset"])
    fixed = reader.read(offset, 30)
    fields = struct.unpack("<4s5H3I2H", fixed)
    signature, _version, flags, compression, _mtime, _mdate, _crc, _compressed, _uncompressed, name_length, extra_length = fields
    require(signature == LOCAL_SIGNATURE, "local ZIP member signature differs")
    require(not flags & 0x0001, "encrypted ZIP member is not admitted")
    require(compression == member["compression"] == 8, "ZIP member compression differs")
    variable = reader.read(offset + 30, name_length + extra_length)
    encoding = "utf-8" if flags & 0x0800 else "cp437"
    try:
        observed_name = variable[:name_length].decode(encoding)
    except UnicodeDecodeError as error:
        raise ScoobyNativeContextReleaseError("local ZIP member name encoding differs") from error
    require(observed_name == member["path"], "local ZIP member name differs")
    payload_offset = offset + 30 + name_length + extra_length
    compressed = reader.read(payload_offset, int(member["compressed_size"]))
    try:
        payload = zlib.decompress(compressed, -zlib.MAX_WBITS)
    except zlib.error as error:
        raise ScoobyNativeContextReleaseError("ZIP member deflate payload differs") from error
    require(len(payload) == member["uncompressed_size"], "ZIP member uncompressed size differs")
    observed_crc = f"{binascii.crc32(payload) & 0xffffffff:08x}"
    require(observed_crc == member["crc32"], "ZIP member CRC32 differs")
    return payload, {
        "path": member["path"],
        "destination": member["destination"],
        "local_header_offset": offset,
        "compression": compression,
        "compressed_size": len(compressed),
        "uncompressed_size": len(payload),
        "crc32": observed_crc,
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def audit(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    config = load_json(config_path.resolve(strict=True), label="native context config")
    require(config.get("schema_version") == SCHEMA, "native context config identity differs")
    require(not output.exists() and not output.is_symlink(), "refusing to overwrite native context output")
    validate_boundary(config)
    artifact_binding = config["authority"]["inventory_artifact"]
    artifact = root / artifact_binding["path"]
    try:
        verify_frozen_tree(artifact)
    except ArtifactError as error:
        raise ScoobyNativeContextReleaseError(f"archive inventory artifact differs: {error}") from error
    require(digest(artifact / "ARTIFACTS.json") == artifact_binding["artifacts_sha256"], "archive inventory identity differs")
    inventory_binding = config["authority"]["inventory"]
    inventory_path = root / inventory_binding["path"]
    require(digest(inventory_path) == inventory_binding["sha256"], "archive inventory file differs")
    inventory = load_json(inventory_path, label="training archive inventory")
    epicardioid_members = validate_inventory(config, inventory)

    reader = HTTPRangeReader(config["archive"]["url"], config["archive"]["size_bytes"])
    output.mkdir(parents=True, mode=0o750)
    records = []
    for member in config["members"]:
        payload, record = extract_member(reader, member)
        destination = Path(member["destination"])
        require(len(destination.parts) == 1 and destination.name == member["destination"], "unsafe output member name")
        write_bytes_exclusive(output / destination, payload, mode=0o640)
        records.append(record)
    receipt = {
        "schema_version": "masld-bench-scooby-epicardioids-native-context-release-receipt-v1",
        "status": "pass_native_context_rows_released_query_encoder_still_absent",
        "config_sha256": digest(config_path),
        "archive": {
            **config["archive"],
            "http_range_request_count": reader.request_count,
            "http_bytes_downloaded": reader.bytes_downloaded,
        },
        "members": records,
        "epicardioid_member_count": len(epicardioid_members),
        "context_encoder_weight_candidates": [],
        "release_boundary": config["release_boundary"],
        "native_context_rows_downloaded": True,
        "native_context_rows_deserialized": False,
        "full_archive_downloaded": False,
        "checkpoint_loaded": False,
        "model_forward_executed": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "terminal_blocker": "released_embedding_rows_exist_but_no_exact_scglue_encoder_guidance_graph_or_query_transform_is_released",
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
