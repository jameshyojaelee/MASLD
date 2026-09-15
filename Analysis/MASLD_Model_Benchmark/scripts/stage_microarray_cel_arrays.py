#!/usr/bin/env python3
"""Stage CEL arrays from a frozen RAW tar onto node-local scratch.

Platform-agnostic successor to ``stage_gse49541_cel_arrays.py``: the series and
its expected record count are arguments rather than constants, so GPL570 and
GPL16686 share one audited code path.

Every decompressed array is checked against the frozen manifest's byte length
and SHA-256 before it is usable, so a summarization job can never read an array
the raw-source audit did not admit.

Staging is label-blind and order-independent: arrays are named by GSM accession
only.  Raw CEL files are never written into a published output-file tree, which
keeps the source's redistribution prohibition intact.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import tarfile

MANIFEST_COLUMNS = (
    "series",
    "sample_accession",
    "member_name",
    "compressed_bytes",
    "decompressed_bytes",
    "decompressed_sha256",
)


class CELStagingError(RuntimeError):
    """Raised when a staged array differs from the frozen raw-source audit."""


def read_manifest(path: Path, *, series: str, expected_records: int) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(MANIFEST_COLUMNS) <= set(reader.fieldnames):
            raise CELStagingError("CEL manifest columns differ")
        rows = [row for row in reader if row["series"] == series]
    if len(rows) != expected_records:
        raise CELStagingError(f"{series} manifest is not {expected_records} records")
    accessions = [row["sample_accession"] for row in rows]
    if len(set(accessions)) != len(accessions):
        raise CELStagingError("CEL manifest repeats a sample accession")
    return sorted(rows, key=lambda row: row["sample_accession"])


def stage_arrays(
    *,
    tar_path: Path,
    manifest_path: Path,
    output: Path,
    series: str,
    expected_records: int,
    platform_id: str,
    accessions: list[str] | None,
) -> dict[str, object]:
    rows = read_manifest(manifest_path, series=series, expected_records=expected_records)
    by_accession = {row["sample_accession"]: row for row in rows}
    requested = sorted(accessions) if accessions else [row["sample_accession"] for row in rows]
    unknown = sorted(set(requested) - set(by_accession))
    if unknown:
        raise CELStagingError(f"requested array is not in the frozen manifest: {unknown}")

    output.mkdir(parents=True, exist_ok=True)
    staged: list[dict[str, object]] = []
    with tarfile.open(tar_path, "r:") as archive:
        members = {member.name: member for member in archive.getmembers()}
        for accession in requested:
            row = by_accession[accession]
            member = members.get(row["member_name"])
            if member is None:
                raise CELStagingError(f"tar member is absent: {row['member_name']}")
            if member.size != int(row["compressed_bytes"]):
                raise CELStagingError(f"compressed size differs for {accession}")
            handle = archive.extractfile(member)
            if handle is None:
                raise CELStagingError(f"tar member is not a regular file: {member.name}")
            destination = output / f"{accession}.CEL"
            digest = sha256()
            written = 0
            with gzip.GzipFile(fileobj=handle) as source, destination.open("wb") as sink:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                    digest.update(block)
                    sink.write(block)
                    written += len(block)
            if written != int(row["decompressed_bytes"]):
                raise CELStagingError(f"decompressed size differs for {accession}")
            observed = digest.hexdigest()
            if observed != row["decompressed_sha256"]:
                raise CELStagingError(f"decompressed SHA-256 differs for {accession}")
            staged.append(
                {
                    "sample_accession": accession,
                    "member_name": row["member_name"],
                    "decompressed_bytes": written,
                    "decompressed_sha256": observed,
                }
            )

    return {
        "schema_version": "masld-bench-microarray-cel-staging-v1",
        "status": "pass_frozen_cel_hashes",
        "series": series,
        "platform_id": platform_id,
        "manifest_records": len(rows),
        "staged_arrays": len(staged),
        "arrays": staged,
        "every_array_matches_frozen_sha256": True,
        "raw_CEL_written_into_published_artifact": False,
        "labels_read": False,
        "GEO_series_matrix_read": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tar", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--series", required=True)
    parser.add_argument("--expected-records", type=int, required=True)
    parser.add_argument("--platform-id", required=True)
    parser.add_argument("--accession", action="append", dest="accessions")
    parser.add_argument("--receipt", type=Path, required=True)
    arguments = parser.parse_args()
    receipt = stage_arrays(
        tar_path=arguments.tar,
        manifest_path=arguments.manifest,
        output=arguments.output,
        series=arguments.series,
        expected_records=arguments.expected_records,
        platform_id=arguments.platform_id,
        accessions=arguments.accessions,
    )
    arguments.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: receipt[key] for key in receipt if key != "arrays"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
