#!/usr/bin/env python3
"""Acquire and integrity-audit raw CEL archives for two transfer cohorts."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import re
import struct
import tarfile
import urllib.request


SOURCES = {
    "GSE49541": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49541/suppl/GSE49541_RAW.tar",
        "bytes": 309_319_680,
        "cel_records": 72,
    },
    "GSE83452": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83452/suppl/GSE83452_RAW.tar",
        "bytes": 1_819_432_960,
        "cel_records": 231,
    },
}
GSM_PATTERN = re.compile(r"(GSM\d+).*\.CEL\.gz$", re.IGNORECASE)


class MicroarrayRawAuditError(RuntimeError):
    """Raised when a raw GEO archive or CEL member does not meet its requirements."""


def identify_cel_format(prefix: bytes) -> str:
    """Identify a CEL container from its documented file signature."""
    if prefix.startswith(b"[CEL]\n") or prefix.startswith(b"[CEL]\r\n"):
        if b"Version=3" not in prefix:
            raise MicroarrayRawAuditError("text CEL lacks version 3 header")
        return "text_v3"
    if len(prefix) >= 8 and struct.unpack("<ii", prefix[:8]) == (64, 4):
        return "xda_v4"
    if len(prefix) >= 2 and prefix[:2] == bytes((59, 1)):
        return "calvin_v1"
    raise MicroarrayRawAuditError("decompressed member lacks a supported CEL signature")


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_accessions(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "sample_accession" not in (reader.fieldnames or ()):
            raise MicroarrayRawAuditError("participant table lacks sample accession")
        values = {row["sample_accession"] for row in reader}
    if not values:
        raise MicroarrayRawAuditError("participant table is empty")
    return values


def download(url: str, path: Path, expected_bytes: int) -> dict[str, object]:
    request = urllib.request.Request(
        url, headers={"User-Agent": "masld-bench-microarray-raw-audit/1.0"}
    )
    digest = sha256()
    observed = 0
    with urllib.request.urlopen(request, timeout=600) as response, path.open("xb") as handle:
        declared = response.headers.get("Content-Length")
        if declared is None or int(declared) != expected_bytes:
            raise MicroarrayRawAuditError("raw archive Content-Length differs")
        while True:
            block = response.read(8 * 1024 * 1024)
            if not block:
                break
            handle.write(block)
            digest.update(block)
            observed += len(block)
        headers = {
            key.lower(): value
            for key in ("Content-Length", "Last-Modified", "ETag")
            if (value := response.headers.get(key)) is not None
        }
    if observed != expected_bytes:
        raise MicroarrayRawAuditError("raw archive byte count differs")
    return {
        "url": url,
        "bytes": observed,
        "sha256": digest.hexdigest(),
        "http_headers": headers,
    }


def audit_tar(path: Path, expected_accessions: set[str]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with tarfile.open(path, mode="r:") as archive:
        members = archive.getmembers()
        if not members:
            raise MicroarrayRawAuditError("raw archive is empty")
        for member in members:
            name = PurePosixPath(member.name)
            if name.is_absolute() or ".." in name.parts or not member.isfile():
                raise MicroarrayRawAuditError("raw archive contains an unsafe member")
            match = GSM_PATTERN.search(name.name)
            if match is None:
                raise MicroarrayRawAuditError(f"non-CEL raw member: {member.name}")
            sample_accession = match.group(1).upper()
            extracted = archive.extractfile(member)
            if extracted is None:
                raise MicroarrayRawAuditError("CEL member cannot be streamed")
            digest = sha256()
            decompressed_bytes = 0
            prefix = b""
            try:
                with gzip.GzipFile(fileobj=extracted, mode="rb") as handle:
                    while True:
                        block = handle.read(4 * 1024 * 1024)
                        if not block:
                            break
                        if not prefix:
                            prefix = block[:32]
                        digest.update(block)
                        decompressed_bytes += len(block)
            except (OSError, EOFError) as error:
                raise MicroarrayRawAuditError(
                    f"CEL gzip integrity failed: {member.name}"
                ) from error
            if decompressed_bytes == 0 or not prefix:
                raise MicroarrayRawAuditError("decompressed CEL is empty")
            cel_format = identify_cel_format(prefix)
            records.append(
                {
                    "sample_accession": sample_accession,
                    "member_name": member.name,
                    "compressed_bytes": member.size,
                    "decompressed_bytes": decompressed_bytes,
                    "decompressed_sha256": digest.hexdigest(),
                    "header_hex_32": prefix.hex(),
                    "cel_format": cel_format,
                }
            )
    observed_accessions = [str(record["sample_accession"]) for record in records]
    if len(set(observed_accessions)) != len(records):
        raise MicroarrayRawAuditError("raw archive repeats one GSM CEL member")
    if set(observed_accessions) != expected_accessions:
        raise MicroarrayRawAuditError("raw CEL and deposited sample axes differ")
    return sorted(records, key=lambda record: str(record["sample_accession"]))


def write_records(path: Path, series: str, rows: list[dict[str, object]]) -> None:
    fields = ("series", *tuple(rows[0]))
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({"series": series, **row})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--participants", type=Path, required=True)
    parser.add_argument("--participants-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if (
        sha256_file(arguments.participants / "ARTIFACTS.json")
        != arguments.participants_artifacts_sha256
    ):
        raise MicroarrayRawAuditError("participant ARTIFACTS SHA-256 differs")
    manifest = json.loads(
        (arguments.participants / "ARTIFACTS.json").read_text(encoding="utf-8")
    )
    if (
        manifest["metadata"].get("participant_join_authoritative") is not True
        or manifest["metadata"].get("model_training_activated") is not False
    ):
        raise MicroarrayRawAuditError("participant artifact role differs")
    accessions = {
        "GSE49541": read_accessions(arguments.participants / "gse49541_participants.tsv"),
        "GSE83452": read_accessions(arguments.participants / "gse83452_records.tsv"),
    }
    arguments.output.mkdir(parents=True, exist_ok=False)
    raw = arguments.output / "raw"
    raw.mkdir()
    receipts: dict[str, object] = {}
    summaries: dict[str, object] = {}
    for series, contract in SOURCES.items():
        archive_path = raw / f"{series}_RAW.tar"
        receipt = download(str(contract["url"]), archive_path, int(contract["bytes"]))
        records = audit_tar(archive_path, accessions[series])
        if len(records) != int(contract["cel_records"]):
            raise MicroarrayRawAuditError(f"{series} CEL record count differs")
        write_records(arguments.output / f"{series}_cel_manifest.tsv", series, records)
        receipts[series] = receipt
        summaries[series] = {
            "cel_records": len(records),
            "compressed_bytes": sum(int(row["compressed_bytes"]) for row in records),
            "decompressed_bytes": sum(int(row["decompressed_bytes"]) for row in records),
            "all_member_gzip_streams_pass": True,
            "cel_formats": dict(
                sorted(Counter(str(row["cel_format"]) for row in records).items())
            ),
            "sample_axis_exact": True,
        }
    summary = {
        "schema_version": "masld-bench-microarray-transfer-raw-audit-v2",
        "status": "pass_raw_CEL_integrity",
        "cohorts": summaries,
        "source_receipts": receipts,
        "raw_array_acquisition_complete": True,
        "raw_CEL_integrity_complete": True,
        "sample_axis_exact": True,
        "platform_probe_mapping_complete": False,
        "normalization_run": False,
        "rights_audit_complete": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }
    (arguments.output / "raw_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
