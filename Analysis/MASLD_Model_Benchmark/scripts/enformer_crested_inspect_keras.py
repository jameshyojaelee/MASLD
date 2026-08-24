#!/usr/bin/env python3
"""Safely inspect the admitted CREsted Enformer Keras archive without loading it."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import io
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import tarfile
import tempfile
from typing import Any, BinaryIO
import zipfile


OUTER_SHA256 = "628f67f540304d4d0e143176dc824ed72b3413f78f2fa2efe5d4f0ab51ea1bcc"
KERAS_NAME = "enformer_crested_human.keras"
KERAS_SHA256 = "29dc3835c1d13a6c6bd93315b554075107cc409a9179b11ee0e00180f67cef56"
KERAS_SIZE = 985_569_856
CLASSES_NAME = "enformer_human_output_classes.tsv"
CLASSES_SHA256 = "83c67554af8900075ee542b1cbf247bf9905fc1b1bc59d5cff11ce64ef58dffb"
CLASSES_SIZE = 230_186


class KerasInspectionError(ValueError):
    """Raised when the restricted converted checkpoint violates its safe contract."""


def _safe_name(value: str) -> str:
    if not value or "\x00" in value or "\\" in value or value.startswith("/"):
        raise KerasInspectionError("archive member has an unsafe name")
    parts = PurePosixPath(value).parts
    if not parts or any(part in {"", ".", ".."} for part in parts):
        raise KerasInspectionError("archive member has unsafe path traversal")
    return PurePosixPath(*parts).as_posix()


def _hash_file(path: Path) -> tuple[str, int]:
    digest = sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def _copy_and_hash(source: BinaryIO, destination: BinaryIO) -> tuple[str, int]:
    digest = sha256()
    size = 0
    for block in iter(lambda: source.read(1024 * 1024), b""):
        destination.write(block)
        digest.update(block)
        size += len(block)
    return digest.hexdigest(), size


def _collect_serialized_classes(value: Any, rows: set[tuple[str, str, str]]) -> None:
    if isinstance(value, dict):
        if "class_name" in value:
            rows.add(
                (
                    str(value.get("module", "")),
                    str(value.get("class_name", "")),
                    str(value.get("registered_name", "")),
                )
            )
        for child in value.values():
            _collect_serialized_classes(child, rows)
    elif isinstance(value, list):
        for child in value:
            _collect_serialized_classes(child, rows)


def _inspect_keras(path: Path) -> tuple[dict[str, Any], dict[str, Any], bytes, bytes]:
    expected_members = {"config.json", "metadata.json", "model.weights.h5"}
    records: list[dict[str, Any]] = []
    captured: dict[str, bytes] = {}
    seen: set[str] = set()
    uncompressed = 0
    with zipfile.ZipFile(path, mode="r") as handle:
        for member in handle.infolist():
            name = _safe_name(member.filename)
            if name in seen:
                raise KerasInspectionError("Keras archive has duplicate members")
            seen.add(name)
            if len(seen) > 32:
                raise KerasInspectionError("Keras archive member ceiling exceeded")
            mode = member.external_attr >> 16
            if stat.S_ISLNK(mode):
                raise KerasInspectionError("Keras archive has a symbolic link")
            if member.is_dir():
                raise KerasInspectionError("Keras archive has an unexpected directory")
            uncompressed += member.file_size
            if uncompressed > 1_500_000_000:
                raise KerasInspectionError("Keras archive expanded-size ceiling exceeded")
            digest = sha256()
            observed = 0
            chunks: list[bytes] = []
            with handle.open(member, mode="r") as source:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(block)
                    observed += len(block)
                    if name in {"config.json", "metadata.json"}:
                        chunks.append(block)
            if observed != member.file_size:
                raise KerasInspectionError("Keras member size differs while streaming")
            if chunks:
                captured[name] = b"".join(chunks)
            records.append(
                {
                    "path": name,
                    "sha256": digest.hexdigest(),
                    "size_bytes": observed,
                    "compressed_size_bytes": member.compress_size,
                }
            )
    if seen != expected_members:
        raise KerasInspectionError("Keras archive member set differs")
    try:
        config = json.loads(captured["config.json"])
        metadata = json.loads(captured["metadata.json"])
    except (KeyError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise KerasInspectionError("Keras JSON metadata is invalid") from error
    classes: set[tuple[str, str, str]] = set()
    _collect_serialized_classes(config, classes)
    dangerous = sorted(
        row for row in classes if row[1].lower() in {"lambda", "pythonlambda"}
    )
    if dangerous:
        raise KerasInspectionError("Keras config contains a Lambda layer")
    inventory = {
        "schema_version": "masld-bench-safe-keras-inventory-v1",
        "status": "pass",
        "member_count": len(records),
        "uncompressed_bytes": uncompressed,
        "archive_opened_by_zipfile_only": True,
        "keras_deserialization_performed": False,
        "members": sorted(records, key=lambda row: row["path"]),
        "serialized_classes": [
            {"module": module, "class_name": name, "registered_name": registered}
            for module, name, registered in sorted(classes)
        ],
        "lambda_layer_present": False,
    }
    return inventory, metadata, captured["config.json"], captured["metadata.json"]


def inspect(archive: Path, output: Path, scratch: Path | None = None) -> dict[str, Any]:
    if output.exists() or archive.is_symlink() or not archive.is_file():
        raise KerasInspectionError("input or output contract is invalid")
    digest, size = _hash_file(archive)
    if digest != OUTER_SHA256 or size != 915_174_055:
        raise KerasInspectionError("outer archive identity differs")
    output.mkdir(parents=True, mode=0o750)
    temporary_parent = None if scratch is None else str(scratch.resolve(strict=True))
    classes_payload: bytes | None = None
    keras_inventory: dict[str, Any] | None = None
    config_payload: bytes | None = None
    metadata_payload: bytes | None = None
    observed_names: set[str] = set()
    with tempfile.TemporaryDirectory(prefix="enformer-crested-", dir=temporary_parent) as work:
        keras_path = Path(work) / KERAS_NAME
        with tarfile.open(archive, mode="r|gz") as handle:
            for member in handle:
                name = _safe_name(member.name)
                if name in observed_names:
                    raise KerasInspectionError("outer archive has duplicate members")
                observed_names.add(name)
                if not member.isfile() or member.issym() or member.islnk():
                    raise KerasInspectionError("outer archive member type differs")
                source = handle.extractfile(member)
                if source is None:
                    raise KerasInspectionError("outer archive member is unreadable")
                if name == KERAS_NAME:
                    with keras_path.open("xb") as destination:
                        member_digest, member_size = _copy_and_hash(source, destination)
                    if member_digest != KERAS_SHA256 or member_size != KERAS_SIZE:
                        raise KerasInspectionError("Keras member identity differs")
                elif name == CLASSES_NAME:
                    buffer = io.BytesIO()
                    member_digest, member_size = _copy_and_hash(source, buffer)
                    if member_digest != CLASSES_SHA256 or member_size != CLASSES_SIZE:
                        raise KerasInspectionError("output-class member identity differs")
                    classes_payload = buffer.getvalue()
                else:
                    raise KerasInspectionError("outer archive member set differs")
        if observed_names != {KERAS_NAME, CLASSES_NAME} or classes_payload is None:
            raise KerasInspectionError("outer archive member set differs")
        keras_inventory, metadata, config_payload, metadata_payload = _inspect_keras(
            keras_path
        )

    try:
        text = classes_payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise KerasInspectionError("output class table is not UTF-8") from error
    rows = text.splitlines()
    if len(rows) != 5313 or any(not row or "\x00" in row for row in rows):
        raise KerasInspectionError(
            "output class table geometry differs: "
            f"rows={len(rows)}"
        )
    label_counts = Counter(rows)
    ordered_label_sha256 = sha256(
        ("\n".join(rows) + "\n").encode("utf-8")
    ).hexdigest()

    assert keras_inventory is not None and config_payload is not None
    assert metadata_payload is not None
    (output / "keras_inventory.json").write_text(
        json.dumps(keras_inventory, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    (output / "keras_config.json").write_bytes(config_payload)
    (output / "keras_metadata.json").write_bytes(metadata_payload)
    (output / "output_classes.txt").write_text(text, encoding="utf-8")
    receipt = {
        "schema_version": "masld-bench-enformer-crested-safe-keras-inspection-v1",
        "status": "pass",
        "model_id": "enformer_crested_restricted_port",
        "outer_archive_sha256": digest,
        "keras_member_sha256": KERAS_SHA256,
        "output_class_table_sha256": CLASSES_SHA256,
        "output_class_rows": len(rows),
        "output_class_format": "headerless_one_UTF8_label_per_line",
        "output_class_unique_labels": len(label_counts),
        "output_class_duplicated_label_values": sum(
            count > 1 for count in label_counts.values()
        ),
        "output_class_maximum_label_multiplicity": max(label_counts.values()),
        "ordered_output_label_sha256": ordered_label_sha256,
        "description_only_track_selection_allowed": False,
        "canonical_5313_track_positional_crosswalk_required": True,
        "safe_outer_member_copy_to_ephemeral_scratch": True,
        "ephemeral_checkpoint_deleted": True,
        "keras_deserialization_performed": False,
        "model_forward_executed": False,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_features_or_labels_read": False,
        "keras_version_from_metadata": metadata.get("keras_version"),
        "serialized_class_count": len(keras_inventory["serialized_classes"]),
        "lambda_layer_present": False,
        "native_sonnet_identity": False,
        "open_champion_eligible": False,
        "allowed_next_action": "canonical_positional_track_crosswalk_and_isolated_safe_mode_Keras_load_with_frozen_custom_object_allowlist",
        "terminal_disposition": "safe_nested_inventory_passed_positional_track_crosswalk_deserialization_and_forward_pending",
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path)
    arguments = parser.parse_args()
    result = inspect(arguments.archive, arguments.output, arguments.scratch)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
