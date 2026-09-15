#!/usr/bin/env python3
"""Resume and verify one participant-contained GSE105127 RNA bundle."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import time
import urllib.error
import urllib.request

from masld_bench.artifacts import freeze_tree, publish_directory_noreplace, verify_frozen_tree
from scripts.download_gse105127_rna_bundle import (
    GSE105127DownloadError,
    audit_fastq,
    digest,
    read_rows,
)


REMOTE_IDENTITY_ARTIFACTS_SHA256 = "7b231c245d2253edd2f7f46fe7ab7040d84d989eb93ccc1c8e9302cd200c0b38"
CONTENT_RANGE = re.compile(r"^bytes (?P<start>[0-9]+)-(?P<end>[0-9]+)/(?P<total>[0-9]+)$")


def validate_remote_identity(root: Path) -> None:
    verify_frozen_tree(root)
    if digest(root / "ARTIFACTS.json", "sha256") != REMOTE_IDENTITY_ARTIFACTS_SHA256:
        raise GSE105127DownloadError("RNA remote identity ARTIFACTS SHA-256 differs")
    receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
    expected = {
        "planned_rows": 57,
        "remote_head_complete_rows": 57,
        "remote_bytes_match_rows": 57,
        "remote_bytes_differ_rows": 0,
        "remote_total_unavailable_rows": 0,
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
    }
    if (
        receipt.get("schema_version") != "masld-bench-gse105127-rna-remote-identity-census-v1"
        or receipt.get("status") != "complete_outcome_free_remote_identity_census"
        or any(receipt.get(key) != value for key, value in expected.items())
    ):
        raise GSE105127DownloadError("RNA remote identity receipt differs")


def validate_response_headers(response: object, *, offset: int, expected_bytes: int) -> None:
    status = int(getattr(response, "status"))
    headers = getattr(response, "headers")
    declared = headers.get("Content-Length")
    if offset == 0:
        if status != 200 or (declared is not None and int(declared) != expected_bytes):
            raise GSE105127DownloadError("initial RNA response identity differs")
        return
    match = CONTENT_RANGE.fullmatch(headers.get("Content-Range", ""))
    if status != 206 or match is None:
        raise GSE105127DownloadError("resumed RNA response is not an exact byte range")
    start, end, total = (int(match.group(key)) for key in ("start", "end", "total"))
    if start != offset or end != expected_bytes - 1 or total != expected_bytes:
        raise GSE105127DownloadError("resumed RNA content range differs")
    if declared is not None and int(declared) != expected_bytes - offset:
        raise GSE105127DownloadError("resumed RNA Content-Length differs")


def resumable_download(url: str, temporary: Path, expected_bytes: int) -> dict[str, object]:
    errors = []
    resumed_requests = 0
    etag = None
    for attempt in range(1, 21):
        offset = temporary.stat().st_size if temporary.exists() else 0
        if offset > expected_bytes:
            raise GSE105127DownloadError("RNA partial file exceeds declared bytes")
        if offset == expected_bytes:
            return {
                "attempts": attempt - 1,
                "resumed_requests": resumed_requests,
                "etag": etag,
            }
        headers = {"User-Agent": "masld-bench-rna-resumable/1.0"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
            resumed_requests += 1
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                validate_response_headers(response, offset=offset, expected_bytes=expected_bytes)
                observed_etag = response.headers.get("ETag")
                if etag is not None and observed_etag is not None and observed_etag != etag:
                    raise GSE105127DownloadError("RNA remote ETag changed during resume")
                etag = etag or observed_etag
                mode = "ab" if offset else "xb"
                with temporary.open(mode) as handle:
                    while True:
                        block = response.read(8 * 1024 * 1024)
                        if not block:
                            break
                        handle.write(block)
                        offset += len(block)
                        if offset > expected_bytes:
                            raise GSE105127DownloadError("RNA FASTQ exceeds declared bytes")
            if offset == expected_bytes:
                return {
                    "attempts": attempt,
                    "resumed_requests": resumed_requests,
                    "etag": etag,
                }
            errors.append(f"short response ended at {offset} of {expected_bytes} bytes")
        except (OSError, urllib.error.URLError, ValueError) as error:
            errors.append(f"{type(error).__name__}: {error}"[:500])
        if attempt < 20:
            time.sleep(min(30, attempt * 2))
    raise GSE105127DownloadError(
        "RNA FASTQ resumable download failed after 20 attempts: " + " | ".join(errors[-3:])
    )


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
    stage.mkdir(parents=True, exist_ok=False)
    path = stage / "reads.fastq.gz"
    temporary = stage / "reads.fastq.gz.part"
    transfer = resumable_download(row["fastq_url"], temporary, int(row["fastq_bytes"]))
    os.replace(temporary, path)
    if digest(path, "md5") != row["fastq_md5"].lower():
        raise GSE105127DownloadError("RNA FASTQ ENA MD5 differs")
    fastq_audit = audit_fastq(
        path,
        expected_reads=int(row["read_count"]),
        expected_bases=int(row["base_count"]),
        read_length=int(row["read_length"]),
    )
    receipt = {
        "schema_version": "masld-bench-gse105127-rna-fastq-v2-resumable",
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
        "remote_identity_artifacts_sha256": REMOTE_IDENTITY_ARTIFACTS_SHA256,
        "transfer": transfer,
        **fastq_audit,
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    (stage / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(stage, {"artifact_class": "gse105127_rna_fastq", "row_id": row["row_id"], "status": "passed"})
    verify_frozen_tree(stage)
    final.parent.mkdir(parents=True, exist_ok=True)
    publish_directory_noreplace(stage, final)
    verify_frozen_tree(final)
    return receipt


def run_bundle(
    *, plan_root: Path, remote_identity_root: Path, bundle_id: int, output: Path, attempt: Path,
) -> dict[str, object]:
    verify_frozen_tree(plan_root)
    validate_remote_identity(remote_identity_root)
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
        (bundle_stage / "receipt.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        freeze_tree(bundle_stage, {"artifact_class": "gse105127_rna_download_bundle", "bundle_id": bundle_id, "status": "passed"})
        publish_directory_noreplace(bundle_stage, bundle)
        verify_frozen_tree(bundle)
    if len(receipts) != len(rows):
        raise GSE105127DownloadError("RNA bundle receipt count differs")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--remote-identity-root", type=Path, required=True)
    parser.add_argument("--bundle-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attempt", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_bundle(
        plan_root=args.plan_root,
        remote_identity_root=args.remote_identity_root,
        bundle_id=args.bundle_id,
        output=args.output,
        attempt=args.attempt,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
