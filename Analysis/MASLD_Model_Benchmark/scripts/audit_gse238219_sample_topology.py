#!/usr/bin/env python3
"""Inventory one GSE238219 replicate without reading expression outcomes."""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path, PurePosixPath
import tarfile
from typing import BinaryIO, Any


SAMPLES = {
    "GSM7660623": "SAMN36706212",
    "GSM7660624": "SAMN36706211",
    "GSM7660625": "SAMN36706210",
    "GSM7660626": "SAMN36706209",
    "GSM7660627": "SAMN36706208",
}


class TopologyAuditError(RuntimeError):
    """Raised when source archive or sample topology differs."""


def file_digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def safe_name(value: str) -> None:
    pure = PurePosixPath(value)
    if pure.is_absolute() or ".." in pure.parts:
        raise TopologyAuditError(f"unsafe archive member: {value}")


def count_text(stream: BinaryIO, *, gzip_encoded: bool) -> dict[str, Any]:
    reader: BinaryIO = gzip.GzipFile(fileobj=stream) if gzip_encoded else stream
    header: str | None = None
    rows = 0
    matrix_shape: list[int] | None = None
    for raw in reader:
        line = raw.decode("utf-8", errors="strict").rstrip("\r\n")
        if not line or line.startswith("%"):
            continue
        if header is None:
            header = line
            if len(line.split()) == 3 and all(value.isdigit() for value in line.split()):
                matrix_shape = [int(value) for value in line.split()]
                break
            continue
        rows += 1
    return {"first_data_line_or_header": header, "data_rows_after_header": rows, "matrix_shape": matrix_shape}


def nested_inventory(stream: BinaryIO) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with tarfile.open(fileobj=stream, mode="r:gz") as archive:
        for member in archive:
            safe_name(member.name)
            record: dict[str, Any] = {
                "name": member.name,
                "size_bytes": member.size,
                "type": "file" if member.isfile() else "directory" if member.isdir() else "other",
            }
            lower = member.name.lower()
            inspect_text = member.isfile() and lower.endswith((".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz", ".mtx.gz"))
            if inspect_text:
                extracted = archive.extractfile(member)
                if extracted is None:
                    raise TopologyAuditError("cannot read nested text member")
                record["text_schema"] = count_text(extracted, gzip_encoded=lower.endswith(".gz"))
            records.append(record)
    return records


def audit(source: Path, sample: str, output: Path) -> dict[str, Any]:
    if sample not in SAMPLES or output.exists() or not source.is_file():
        raise TopologyAuditError("sample topology request differs")
    output.mkdir(parents=True, exist_ok=False)
    selected: list[dict[str, Any]] = []
    with tarfile.open(source, mode="r:") as archive:
        members = [member for member in archive if member.isfile() and member.name.startswith(sample)]
        if len(members) != 2:
            raise TopologyAuditError(f"expected two sample archives, observed {len(members)}")
        for member in sorted(members, key=lambda value: value.name):
            safe_name(member.name)
            extracted = archive.extractfile(member)
            if extracted is None:
                raise TopologyAuditError("cannot read selected sample archive")
            selected.append(
                {
                    "name": member.name,
                    "size_bytes": member.size,
                    "members": nested_inventory(extracted),
                }
            )
    result = {
        "schema_version": "masld-bench-gse238219-sample-topology-audit-v1",
        "sample": sample,
        "biosample": SAMPLES[sample],
        "source_path": source.as_posix(),
        "source_sha256": file_digest(source),
        "selected_archives": selected,
        "expression_values_read": False,
        "matrix_dimensions_read": True,
        "guide_assignment_schema_read": True,
        "biological_unit_status": "deposited_replicate_label_not_yet_independence_proof",
        "gse313774_accessed": False,
    }
    (output / "sample_topology.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"sample": sample, "archives": len(selected)}, sort_keys=True))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--sample", choices=tuple(SAMPLES), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.source, args.sample, args.output)


if __name__ == "__main__":
    main()
