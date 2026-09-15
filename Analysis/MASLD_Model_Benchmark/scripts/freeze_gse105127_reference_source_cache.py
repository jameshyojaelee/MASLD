#!/usr/bin/env python3
"""Freeze exact GSE105127 reference source bytes from a preserved attempt."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from scripts.freeze_gse105127_reference_bundle import (
    FORWARD_MD5,
    GSE105127ReferenceError,
    REVERSE_MD5,
    SOURCE_BYTES,
    digest_file,
)
from scripts.freeze_gse105127_reference_bundle_legacy_exact import SOURCE_SHA256_COMPRESSED


def build(*, source_root: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise GSE105127ReferenceError(f"refusing to overwrite source cache: {output}")
    expected = {
        "human_g1k_v37.fasta.gz": ("sha256", SOURCE_SHA256_COMPRESSED, SOURCE_BYTES),
        "hg19ToHg38.over.chain.gz": ("md5", FORWARD_MD5, 227_698),
        "hg38ToHg19.over.chain.gz": ("md5", REVERSE_MD5, 1_246_411),
    }
    for name, (algorithm, digest, size) in expected.items():
        source = source_root / name
        if not source.is_file() or source.stat().st_size != size:
            raise GSE105127ReferenceError(f"reference source cache input size differs: {name}")
        if digest_file(source, algorithm) != digest:
            raise GSE105127ReferenceError(f"reference source cache input digest differs: {name}")
    for name in ("hg19ToHg38.over.chain.gz", "hg38ToHg19.over.chain.gz"):
        if subprocess.run(["gzip", "-t", str(source_root / name)], check=False).returncode != 0:
            raise GSE105127ReferenceError(f"reference source cache gzip differs: {name}")

    output.mkdir(parents=True)
    rows = []
    for name, (algorithm, digest, size) in expected.items():
        target = output / name
        with (source_root / name).open("rb") as source, target.open("xb") as destination:
            shutil.copyfileobj(source, destination, length=8 * 1024 * 1024)
        if target.stat().st_size != size or digest_file(target, algorithm) != digest:
            raise GSE105127ReferenceError(f"copied reference source differs: {name}")
        rows.append({"name": name, "bytes": size, "algorithm": algorithm, "digest": digest})
    receipt = {
        "schema_version": "masld-bench-gse105127-reference-source-cache-v1",
        "status": "passed_exact_source_bytes",
        "files": rows,
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(output, {
        "artifact_class": "gse105127_reference_source_cache",
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "status": "passed",
    })
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(source_root=args.source_root, output=args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
