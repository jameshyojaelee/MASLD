#!/usr/bin/env python3
"""Download and verify one participant-contained GSE105127 RNA bundle."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import md5, sha256
import json
import os
from pathlib import Path
import time
import urllib.request

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace, verify_frozen_tree


class GSE105127DownloadError(RuntimeError):
    """Raised when a raw RNA member differs from ENA."""


def read_rows(path: Path, bundle_id: int) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        forbidden = {"phenotype", "label", "outcome", "disease", "fibrosis", "nas", "sex", "age", "bmi", "outer_fold"}
        if forbidden & fields:
            raise GSE105127DownloadError("RNA plan violates the label firewall")
        rows = [row for row in reader if int(row["bundle_id"]) == bundle_id]
    if not rows:
        raise GSE105127DownloadError("RNA bundle is empty")
    groups = {row["participant_group_id"] for row in rows}
    if any(sum(member["participant_group_id"] == group for member in rows) != 3 for group in groups):
        raise GSE105127DownloadError("participant zones are split inside RNA bundle")
    return rows


def digest(path: Path, algorithm: str) -> str:
    value = {"md5": md5, "sha256": sha256}[algorithm]()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def audit_fastq(path: Path, expected_reads: int, expected_bases: int, read_length: int) -> dict[str, int]:
    reads = 0
    bases = 0
    try:
        with gzip.open(path, "rb") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().rstrip(b"\r\n")
                plus = handle.readline()
                quality = handle.readline().rstrip(b"\r\n")
                if (
                    not header.startswith(b"@")
                    or not plus.startswith(b"+")
                    or len(sequence) != read_length
                    or len(quality) != len(sequence)
                ):
                    raise GSE105127DownloadError("RNA FASTQ record structure differs")
                reads += 1
                bases += len(sequence)
    except (EOFError, OSError) as error:
        raise GSE105127DownloadError("RNA FASTQ gzip stream differs") from error
    if reads != expected_reads or bases != expected_bases:
        raise GSE105127DownloadError("RNA FASTQ read/base census differs")
    return {"read_count": reads, "base_count": bases, "read_length": read_length}


def download_row(row: dict[str, str], output: Path, attempt: Path) -> dict[str, object]:
    final = output / "participants" / row["row_id"]
    if (final / "COMPLETE").is_file():
        verify_frozen_tree(final)
        receipt = json.loads((final / "receipt.json").read_text(encoding="utf-8"))
        expected = {
            "row_id": row["row_id"],
            "participant_group_id": row["participant_group_id"],
            "zone": row["zone"],
            "md5": row["fastq_md5"].lower(),
            "read_count": int(row["read_count"]),
            "base_count": int(row["base_count"]),
            "read_length": int(row["read_length"]),
        }
        if any(receipt.get(key) != value for key, value in expected.items()):
            raise GSE105127DownloadError("replayed RNA member differs from plan")
        return receipt
    if final.exists():
        raise GSE105127DownloadError("incomplete final RNA member exists")
    stage = attempt / f"{row['row_id']}.staging"
    if stage.exists():
        raise GSE105127DownloadError("attempt-local RNA stage exists")
    stage.mkdir(parents=True)
    path = stage / "reads.fastq.gz"
    temporary = stage / "reads.fastq.gz.part"
    last_error: Exception | None = None
    for retry in range(1, 4):
        try:
            request = urllib.request.Request(row["fastq_url"], headers={"User-Agent": "masld-bench-rna/1.0"})
            size = 0
            with urllib.request.urlopen(request, timeout=300) as response, temporary.open("xb") as handle:
                while True:
                    block = response.read(8 * 1024 * 1024)
                    if not block:
                        break
                    size += len(block)
                    if size > int(row["fastq_bytes"]):
                        raise GSE105127DownloadError("RNA FASTQ exceeds declared bytes")
                    handle.write(block)
            if size != int(row["fastq_bytes"]):
                raise GSE105127DownloadError("RNA FASTQ bytes differ")
            os.replace(temporary, path)
            last_error = None
            break
        except Exception as error:
            last_error = error
            temporary.unlink(missing_ok=True)
            path.unlink(missing_ok=True)
            if retry < 3:
                time.sleep(retry * 2)
    if last_error is not None:
        raise GSE105127DownloadError("RNA FASTQ download failed after retries") from last_error
    if digest(path, "md5") != row["fastq_md5"].lower():
        raise GSE105127DownloadError("RNA FASTQ ENA MD5 differs")
    fastq_audit = audit_fastq(
        path,
        expected_reads=int(row["read_count"]),
        expected_bases=int(row["base_count"]),
        read_length=int(row["read_length"]),
    )
    receipt = {
        "schema_version": "masld-bench-gse105127-rna-fastq-v1",
        "status": "complete_verified",
        "row_id": row["row_id"],
        "participant_group_id": row["participant_group_id"],
        "zone": row["zone"],
        "pairing_topology": "adjacent_section",
        "run_accession": row["run_accession"],
        "bytes": path.stat().st_size,
        "md5": row["fastq_md5"].lower(),
        "sha256": digest(path, "sha256"),
        "gzip_integrity": True,
        **fastq_audit,
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    (stage / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(stage, {"artifact_class": "gse105127_rna_fastq", "row_id": row["row_id"], "status": "passed"})
    verify_frozen_tree(stage)
    final.parent.mkdir(parents=True, exist_ok=True)
    publish_directory_noreplace(stage, final)
    verify_frozen_tree(final)
    return receipt


def run_bundle(*, plan_root: Path, bundle_id: int, output: Path, attempt: Path) -> dict[str, object]:
    verify_frozen_tree(plan_root)
    rows = read_rows(plan_root / "rna_rows.tsv", bundle_id)
    output.mkdir(parents=True, exist_ok=True)
    attempt.mkdir(parents=True, exist_ok=True)
    receipts = [download_row(row, output, attempt) for row in rows]
    result = {
        "schema_version": "masld-bench-gse105127-rna-download-bundle-v1",
        "status": "complete",
        "bundle_id": bundle_id,
        "participants": len({row["participant_group_id"] for row in rows}),
        "participant_zone_rows": len(rows),
        "bytes": sum(int(row["fastq_bytes"]) for row in rows),
        "all_md5_and_gzip_verified": True,
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    bundle = output / f"bundle_{bundle_id:02d}"
    if (bundle / "COMPLETE").is_file():
        verify_frozen_tree(bundle)
        if json.loads((bundle / "receipt.json").read_text(encoding="utf-8")) != result:
            raise GSE105127DownloadError("replayed RNA bundle receipt differs")
    elif bundle.exists():
        raise GSE105127DownloadError("incomplete final RNA bundle receipt exists")
    else:
        bundle_stage = attempt / f"bundle_{bundle_id:02d}.staging"
        bundle_stage.mkdir(exist_ok=False)
        (bundle_stage / "receipt.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        freeze_tree(bundle_stage, {"artifact_class": "gse105127_rna_download_bundle", "bundle_id": bundle_id, "status": "passed"})
        publish_directory_noreplace(bundle_stage, bundle)
        verify_frozen_tree(bundle)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--bundle-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attempt", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_bundle(plan_root=args.plan_root, bundle_id=args.bundle_id, output=args.output, attempt=args.attempt), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
