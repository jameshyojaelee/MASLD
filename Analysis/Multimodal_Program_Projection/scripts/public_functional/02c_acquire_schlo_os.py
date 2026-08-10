#!/usr/bin/env python3
"""Acquire the checksum-recorded same-study OS-HLO source amendment."""

from __future__ import annotations

import datetime as dt
import gzip
import json
import os
import shutil
import subprocess

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, md5_file, sha256_file, write_tsv


def main() -> None:
    amendment_path = CANDIDATE_ROOT / "SOURCE_AMENDMENT_01.json"
    amendment = json.loads(amendment_path.read_text(encoding="utf-8"))
    root = CANDIDATE_ROOT / "sources/GSE207889"
    compressed = root / amendment["local_name"]
    partial = compressed.with_suffix(compressed.suffix + ".part")
    subprocess.run(["curl", "--fail", "--location", "--retry", "6", "--continue-at", "-", "--output", str(partial), amendment["url"]], check=True)
    os.replace(partial, compressed)
    h5ad = compressed.with_suffix("")
    h5ad_partial = h5ad.with_suffix(h5ad.suffix + ".part")
    with gzip.open(compressed, "rb") as source, h5ad_partial.open("wb") as target:
        shutil.copyfileobj(source, target, length=16 * 1024 * 1024)
    os.replace(h5ad_partial, h5ad)
    row = {
        "dataset_id": "GSE207889", "file_id": "os_h5ad_source_amendment",
        "source_url": amendment["url"], "local_path": str(compressed.relative_to(PROJECT_ROOT)),
        "retrieved_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "size_bytes": compressed.stat().st_size,
        "md5": md5_file(compressed), "sha256": sha256_file(compressed), "status": "downloaded_source_amendment",
        "required_for_gate": "true", "specification_sha256": amendment["parent_specification_sha256"],
        "source_amendment_sha256": amendment["amendment_sha256"],
    }
    write_tsv(root / "source_amendment_manifest.tsv", [row], list(row))
    print(json.dumps(row, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
