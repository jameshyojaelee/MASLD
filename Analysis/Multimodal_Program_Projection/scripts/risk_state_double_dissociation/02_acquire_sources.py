#!/usr/bin/env python3
"""Acquire checksum-addressable Plan 42 public source products after seal validation."""

from __future__ import annotations

import datetime as dt
import json
import os
import tempfile
import urllib.request
from pathlib import Path

from risk_state_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    md5_file,
    require_validated_seal,
    sha256_file,
    write_tsv,
)


SOURCES = [
    {
        "source_id": "CLCC1_ST1",
        "filename": "41586_2025_10064_MOESM3_ESM.xlsx",
        "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-025-10064-4/MediaObjects/41586_2025_10064_MOESM3_ESM.xlsx",
        "expected_size": 4_215_920,
        "expected_md5": "c78f0e43640261bd77f26652b4fd2f93",
    },
    {
        "source_id": "GSE313544_SOFT",
        "filename": "GSE313544_family.soft.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313544/soft/GSE313544_family.soft.gz",
    },
    {
        "source_id": "GSE313544_COUNTS",
        "filename": "GSE313544_normalized_counts.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313544/suppl/GSE313544_normalized_counts.txt.gz",
    },
    {
        "source_id": "GSE158182_SOFT",
        "filename": "GSE158182_family.soft.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158182/soft/GSE158182_family.soft.gz",
    },
    {
        "source_id": "GSE158182_COUNTS",
        "filename": "GSE158182_RawCounts_Subread_Genes.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158182/suppl/GSE158182_RawCounts_Subread_Genes.txt.gz",
    },
    {
        "source_id": "HMSMA_PORTAL_DOWNLOAD_PAGE",
        "filename": "hmsma_download.html",
        "url": "https://db.genomics.cn/stomics/hmsma/download/",
    },
    {
        "source_id": "HMSMA_PORTAL_ACCESS_JS",
        "filename": "hmsma.popup.js",
        "url": "https://db.genomics.cn/stomics/assets/js/hmsma.popup.js",
    },
    {
        "source_id": "HMSMA_SOURCE_CODE",
        "filename": "Spatial_multiomics_analysis_MASLD-85767ada483db0d326477b7b3d5d8d355426664b.tar.gz",
        "url": "https://github.com/OMIC-coding/Spatial_multiomics_analysis_MASLD/archive/85767ada483db0d326477b7b3d5d8d355426664b.tar.gz",
    },
]


def acquire(source: dict[str, object], destination: Path) -> str:
    expected_size = source.get("expected_size")
    expected_md5 = source.get("expected_md5")
    if destination.is_file():
        size_ok = expected_size is None or destination.stat().st_size == expected_size
        md5_ok = expected_md5 is None or md5_file(destination) == expected_md5
        if size_ok and md5_ok:
            return "reused_verified"
        raise RuntimeError(f"Existing source does not match frozen identity: {destination}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        request = urllib.request.Request(
            str(source["url"]), headers={"User-Agent": "MASLD-Plan42-source-gate/1.0"}
        )
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as handle:
            while chunk := response.read(8 * 1024 * 1024):
                handle.write(chunk)
        if expected_size is not None and temporary.stat().st_size != expected_size:
            raise RuntimeError(
                f"Size mismatch for {source['source_id']}: {temporary.stat().st_size} != {expected_size}"
            )
        if expected_md5 is not None and md5_file(temporary) != expected_md5:
            raise RuntimeError(
                f"MD5 mismatch for {source['source_id']}: {md5_file(temporary)} != {expected_md5}"
            )
        os.replace(temporary, destination)
    except Exception:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise
    return "downloaded_verified"


def main() -> None:
    seal = require_validated_seal()
    source_root = CANDIDATE_ROOT / "source"
    rows: list[dict[str, object]] = []
    for source in SOURCES:
        destination = source_root / str(source["source_id"]) / str(source["filename"])
        status = acquire(source, destination)
        rows.append(
            {
                "source_id": source["source_id"],
                "url": source["url"],
                "local_path": str(destination.relative_to(PROJECT_ROOT)),
                "size_bytes": destination.stat().st_size,
                "md5": md5_file(destination),
                "sha256": sha256_file(destination),
                "retrieved_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "status": status,
                "specification_sha256": seal["specification_sha256"],
            }
        )
    write_tsv(
        CANDIDATE_ROOT / "source_manifest.tsv",
        rows,
        [
            "source_id", "url", "local_path", "size_bytes", "md5", "sha256",
            "retrieved_at_utc", "status", "specification_sha256",
        ],
    )
    print(json.dumps(rows, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

