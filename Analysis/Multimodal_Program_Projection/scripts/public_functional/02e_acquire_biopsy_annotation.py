#!/usr/bin/env python3
"""Acquire and extract the checksum-pinned GPL16686 annotation SQLite file."""

from __future__ import annotations

import datetime as dt
import json
import os
import shutil
import subprocess
import tarfile

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, md5_file, sha256_file, write_tsv


def main() -> None:
    amendment = json.loads((CANDIDATE_ROOT / "SOURCE_AMENDMENT_02.json").read_text(encoding="utf-8"))
    root = CANDIDATE_ROOT / "sources/GSE106737/annotation"
    root.mkdir(parents=True, exist_ok=True)
    archive = root / "hugene20sttranscriptcluster.db_8.8.0.tar.gz"
    partial = archive.with_suffix(archive.suffix + ".part")
    subprocess.run(["curl", "--fail", "--location", "--retry", "6", "--continue-at", "-", "--output", str(partial), amendment["url"]], check=True)
    os.replace(partial, archive)
    if md5_file(archive) != amendment["published_md5"]:
        raise RuntimeError("Bioconductor package MD5 does not match the published PACKAGES index")
    with tarfile.open(archive, "r:gz") as tar:
        members = [member for member in tar.getmembers() if member.name.endswith(".sqlite")]
        if len(members) != 1:
            raise RuntimeError(f"Expected one annotation SQLite file, found {len(members)}")
        extracted = tar.extractfile(members[0])
        if extracted is None:
            raise RuntimeError("Unable to read annotation SQLite member")
        sqlite = root / "hugene20sttranscriptcluster.sqlite"
        sqlite_partial = sqlite.with_suffix(".sqlite.part")
        with sqlite_partial.open("wb") as handle:
            shutil.copyfileobj(extracted, handle)
        os.replace(sqlite_partial, sqlite)
    row = {
        "dataset_id": "GSE106737", "file_id": "GPL16686_annotation_sqlite",
        "source_url": amendment["url"], "local_path": str(sqlite.relative_to(PROJECT_ROOT)),
        "retrieved_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "archive_size_bytes": archive.stat().st_size, "archive_md5": md5_file(archive),
        "archive_sha256": sha256_file(archive), "sqlite_size_bytes": sqlite.stat().st_size,
        "sqlite_sha256": sha256_file(sqlite), "source_amendment_sha256": amendment["amendment_sha256"],
    }
    write_tsv(root / "annotation_source_manifest.tsv", [row], list(row))
    print(json.dumps(row, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
