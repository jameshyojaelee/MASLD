#!/usr/bin/env python3
"""Acquire one public source lane with resume, checksums, and fail-closed provenance."""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import os
import shutil
import subprocess
import tarfile
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from public_functional_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_text,
    md5_file,
    read_tsv,
    require_sealed,
    require_within,
    sha256_file,
    write_tsv,
)


DATASETS = {"GSE200418", "GSE207889", "GSE253380", "GSE106737", "CLCC1"}
CLCC1_MD5 = "fd3fcd19feec6572ff7df90c33f9d906"
CLCC1_SIZE = 4_215_858


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    subprocess.run(
        [
            "curl", "--fail", "--location", "--retry", "6", "--retry-delay", "5",
            "--continue-at", "-", "--output", str(partial), url,
        ],
        check=True,
    )
    os.replace(partial, destination)


def recover_clcc1_from_oa(destination: Path) -> str:
    api = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC11185754"
    with urllib.request.urlopen(api, timeout=60) as response:
        root = ET.fromstring(response.read())
    links = [node.attrib.get("href", "") for node in root.findall(".//link")]
    packages = [link for link in links if link.endswith(".tar.gz")]
    if not packages:
        raise RuntimeError("PMC OA API did not return a tar.gz package")
    package_url = packages[0].replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
    archive = destination.parent / "PMC11185754.tar.gz"
    download(package_url, archive)
    with tarfile.open(archive, "r:gz") as tar:
        candidates = [member for member in tar.getmembers() if Path(member.name).name == "media-1.xlsx"]
        if len(candidates) != 1:
            raise RuntimeError(f"Expected one media-1.xlsx in OA package; found {len(candidates)}")
        member = candidates[0]
        extracted = tar.extractfile(member)
        if extracted is None:
            raise RuntimeError("Unable to extract CLCC1 workbook member")
        partial = destination.with_suffix(destination.suffix + ".part")
        with partial.open("wb") as handle:
            shutil.copyfileobj(extracted, handle)
        os.replace(partial, destination)
    return package_url


def materialize_h5ad(compressed: Path) -> Path:
    destination = compressed.with_suffix("")
    if destination.is_file() and destination.stat().st_size > 0:
        return destination
    partial = destination.with_suffix(destination.suffix + ".part")
    with gzip.open(compressed, "rb") as source, partial.open("wb") as target:
        shutil.copyfileobj(source, target, length=16 * 1024 * 1024)
    os.replace(partial, destination)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    args = parser.parse_args()
    seal = require_sealed()
    source_rows = [
        row for row in read_tsv(CANDIDATE_ROOT / "frozen_spec/source_files.tsv")
        if row["dataset_id"] == args.dataset
    ]
    if not source_rows:
        raise RuntimeError(f"No frozen source rows for {args.dataset}")

    raw_root = require_within(CANDIDATE_ROOT / "sources" / args.dataset)
    raw_root.mkdir(parents=True, exist_ok=True)
    unseal = CANDIDATE_ROOT / "UNSEALED.json"
    if not unseal.exists():
        atomic_write_text(
            unseal,
            json.dumps(
                {
                    "unsealed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "reason": "public source acquisition and schema gates began after outcome-blind seal",
                    "specification_sha256": seal["specification_sha256"],
                },
                indent=2,
                sort_keys=True,
            ) + "\n",
        )

    manifest = []
    errors = []
    for row in source_rows:
        destination = raw_root / row["local_name"]
        used_url = row["url"]
        status = "downloaded"
        try:
            if not destination.is_file() or destination.stat().st_size == 0:
                try:
                    download(used_url, destination)
                except Exception:
                    if args.dataset != "CLCC1":
                        raise
                    used_url = recover_clcc1_from_oa(destination)
            if args.dataset == "CLCC1" and (
                destination.stat().st_size != CLCC1_SIZE or md5_file(destination) != CLCC1_MD5
            ):
                try:
                    used_url = recover_clcc1_from_oa(destination)
                except Exception as exc:
                    status = "checksum_failed"
                    errors.append(f"CLCC1 checksum recovery failed: {exc}")
            if row["extract_gzip"].lower() == "true" and destination.is_file():
                materialize_h5ad(destination)
            manifest.append(
                {
                    "dataset_id": args.dataset,
                    "file_id": row["file_id"],
                    "source_url": used_url,
                    "local_path": str(destination.relative_to(PROJECT_ROOT)),
                    "retrieved_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "size_bytes": destination.stat().st_size if destination.exists() else "",
                    "md5": md5_file(destination) if destination.exists() else "",
                    "sha256": sha256_file(destination) if destination.exists() else "",
                    "status": status,
                    "required_for_gate": row["required_for_gate"],
                    "specification_sha256": seal["specification_sha256"],
                }
            )
        except Exception as exc:
            errors.append(f"{row['file_id']}: {type(exc).__name__}: {exc}")
            manifest.append(
                {
                    "dataset_id": args.dataset,
                    "file_id": row["file_id"],
                    "source_url": used_url,
                    "local_path": str(destination.relative_to(PROJECT_ROOT)),
                    "retrieved_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "size_bytes": destination.stat().st_size if destination.exists() else "",
                    "md5": md5_file(destination) if destination.exists() else "",
                    "sha256": sha256_file(destination) if destination.exists() else "",
                    "status": "download_failed",
                    "required_for_gate": row["required_for_gate"],
                    "specification_sha256": seal["specification_sha256"],
                }
            )

    fields = [
        "dataset_id", "file_id", "source_url", "local_path", "retrieved_at_utc",
        "size_bytes", "md5", "sha256", "status", "required_for_gate", "specification_sha256",
    ]
    write_tsv(raw_root / "source_manifest.tsv", manifest, fields)
    atomic_write_text(raw_root / "acquisition_status.json", json.dumps({"dataset_id": args.dataset, "errors": errors, "n_files": len(manifest)}, indent=2, sort_keys=True) + "\n")
    if errors:
        raise SystemExit("; ".join(errors))
    print(f"ACQUISITION_COMPLETE\t{args.dataset}\t{len(manifest)} files")


if __name__ == "__main__":
    main()
