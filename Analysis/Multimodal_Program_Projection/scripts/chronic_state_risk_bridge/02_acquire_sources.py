#!/usr/bin/env python3
"""Acquire frozen Plan 43 public sources after the outcome-blind seal."""

from __future__ import annotations

import argparse
import datetime as dt
import os
import time
import urllib.request
from pathlib import Path

from bridge_common import (
    CANDIDATE_ROOT,
    CONFIG_ROOT,
    read_tsv,
    require_validated_seal,
    sha256_file,
    write_tsv,
)


GROUP_DATASETS = {
    "human": {"GSE83452", "GSE48452"},
    "crop": {"GSE281160", "ZHU_SUPP", "ZHU_SOURCE"},
    "protein": {"PXD052787"},
}


def acquire(url: str, destination: Path, retries: int = 4) -> tuple[int, str]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "MASLD-Plan43-source-audit/1.0 (public-data research)",
                    "Accept": "*/*",
                },
            )
            with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as handle:
                expected = response.headers.get("Content-Length")
                while True:
                    chunk = response.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
            observed = partial.stat().st_size
            if expected is not None and observed != int(expected):
                raise IOError(f"content-length mismatch: {observed} != {expected}")
            os.replace(partial, destination)
            return observed, sha256_file(destination)
        except Exception:
            if attempt == retries:
                raise
            time.sleep(5 * attempt)
    raise AssertionError("unreachable")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", choices=sorted(GROUP_DATASETS), required=True)
    args = parser.parse_args()
    require_validated_seal()

    rows = [
        row
        for row in read_tsv(CONFIG_ROOT / "source_registry.tsv")
        if row["dataset_id"] in GROUP_DATASETS[args.group]
    ]
    if not rows:
        raise RuntimeError(f"No registered sources for {args.group}")
    output: list[dict[str, object]] = []
    for row in rows:
        destination = CANDIDATE_ROOT / "sources" / row["dataset_id"] / row["local_name"]
        record: dict[str, object] = {
            "group": args.group,
            "dataset_id": row["dataset_id"],
            "prospective_status": row["prospective_status"],
            "url": row["url"],
            "relative_path": str(destination.relative_to(CANDIDATE_ROOT)),
            "retrieved_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "acquired": "false",
            "size_bytes": "",
            "sha256": "",
            "error_type": "",
            "error_message": "",
        }
        try:
            if destination.exists():
                raise FileExistsError(f"refusing to overwrite {destination}")
            size, digest = acquire(row["url"], destination)
            record.update(acquired="true", size_bytes=size, sha256=digest)
        except Exception as error:  # source availability is a gate, not a software crash
            record.update(
                error_type=type(error).__name__,
                error_message=str(error).replace("\t", " ").replace("\n", " ")[:1000],
            )
        output.append(record)

    fields = [
        "group", "dataset_id", "prospective_status", "url", "relative_path",
        "retrieved_at_utc", "acquired", "size_bytes", "sha256", "error_type",
        "error_message",
    ]
    write_tsv(CANDIDATE_ROOT / "acquisition" / f"{args.group}_source_manifest.tsv", output, fields)
    print(f"{args.group}: acquired {sum(row['acquired'] == 'true' for row in output)}/{len(output)} registered objects")


if __name__ == "__main__":
    main()
