#!/usr/bin/env python3
"""Stage GSE49541 CEL arrays from the frozen RAW tar onto node-local scratch.

Every decompressed array is checked against the frozen
``GSE49541_cel_manifest.tsv`` byte length and SHA-256 before it is usable, so a
summarization job can never read an array the raw-source audit did not admit.

Staging is label-blind and order-independent: arrays are named by GSM accession
only.  Raw CEL files are never written into a published output-file tree, which
keeps the source's redistribution prohibition intact.
"""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import tarfile

EXPECTED_SERIES = "GSE49541"
EXPECTED_RECORDS = 72
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


def read_manifest(path: Path) -> list[dict[str, str]]:
    import csv

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(MANIFEST_COLUMNS) <= set(reader.fieldnames):
            raise CELStagingError("CEL manifest columns differ")
        rows = list(reader)
    if len(rows) != EXPECTED_RECORDS:
        raise CELStagingError(f"CEL manifest is not {EXPECTED_RECORDS} records")
    if {row["series"] for row in rows} != {EXPECTED_SERIES}:
        raise CELStagingError("CEL manifest carries a foreign series")
    accessions = [row["sample_accession"] for row in rows]
    if len(set(accessions)) != len(accessions):
        raise CELStagingError("CEL manifest repeats a sample accession")
    return sorted(rows, key=lambda row: row["sample_accession"])


def stage_arrays(
    *, tar_path: Path, manifest_path: Path, output: Path, accessions: list[str] | None
) -> dict[str, object]:
    rows = read_manifest(manifest_path)
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

    receipt = {
        "schema_version": "masld-bench-gse49541-cel-staging-v1",
        "status": "pass_frozen_cel_hashes",
        "series": EXPECTED_SERIES,
        "platform_id": "GPL570",
        "manifest_records": len(rows),
        "staged_arrays": len(staged),
        "arrays": staged,
        "every_array_matches_frozen_sha256": True,
        "raw_CEL_written_into_published_artifact": False,
        "labels_read": False,
        "GEO_series_matrix_read": False,
    }
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tar", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--accession", action="append", dest="accessions")
    parser.add_argument("--receipt", type=Path, required=True)
    arguments = parser.parse_args()
    receipt = stage_arrays(
        tar_path=arguments.tar,
        manifest_path=arguments.manifest,
        output=arguments.output,
        accessions=arguments.accessions,
    )
    arguments.receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: receipt[key] for key in receipt if key != "arrays"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
