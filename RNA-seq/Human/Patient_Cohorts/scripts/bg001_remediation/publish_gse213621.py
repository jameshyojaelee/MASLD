#!/usr/bin/env python3
"""Publish a fully validated GSE213621 staging merge, with count renamed last."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path


FILES = ("gene_counts.txt", "gene_counts.txt.summary", "featureCounts.log", "environment.txt", "validation.json")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    if not (root / ".bg001_candidate_root").is_file():
        raise SystemExit("Missing candidate sentinel")
    parent = root / "counts/GSE213621"
    stage = parent / ".merge_staging"
    if not stage.is_dir() or stage.is_symlink():
        raise SystemExit("Missing or symlinked GSE213621 merge staging directory")
    staged = {name: stage / name for name in FILES}
    for path in staged.values():
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise SystemExit(f"Missing, empty, or symlinked staged artifact: {path}")
    validation = json.loads(staged["validation.json"].read_text())
    expected_hashes = {
        "gene_counts.txt": validation.get("counts_sha256"),
        "gene_counts.txt.summary": validation.get("summary_sha256"),
        "featureCounts.log": validation.get("log_sha256"),
        "environment.txt": validation.get("environment_sha256"),
    }
    if validation.get("status") != "PASS" or validation.get("n_genes") != 86_369 or validation.get("n_samples") != 367:
        raise SystemExit("Staged GSE213621 validation payload is not PASS/86,369/367")
    for name, expected in expected_hashes.items():
        if not expected or sha256(staged[name]) != expected:
            raise SystemExit(f"Staged GSE213621 hash differs: {name}")

    publication = parent / "merge_publication.json"
    marker = parent / "MERGE_COMPLETE"
    final = {name: parent / name for name in FILES}
    for path in (*final.values(), publication, marker):
        if path.exists() or path.is_symlink():
            raise SystemExit(f"Refusing existing GSE213621 publication artifact: {path}")
    publication_payload = {
        "status": "VALIDATED_STAGING",
        "run_id": root.name,
        "publish_order": [
            "environment.txt", "featureCounts.log", "gene_counts.txt.summary",
            "validation.json", "gene_counts.txt",
        ],
        "hashes": {name: sha256(path) for name, path in staged.items()},
    }
    publication_tmp = publication.with_suffix(".json.tmp")
    publication_tmp.write_text(json.dumps(publication_payload, indent=2, sort_keys=True) + "\n")
    os.replace(publication_tmp, publication)
    for name in publication_payload["publish_order"]:
        os.rename(staged[name], final[name])
    for name, expected in publication_payload["hashes"].items():
        if sha256(final[name]) != expected:
            raise SystemExit(f"Post-publication GSE213621 hash differs: {name}")
    marker_tmp = marker.with_suffix(".tmp")
    marker_tmp.write_text(sha256(publication) + "\n")
    os.replace(marker_tmp, marker)
    print("PASS: published fully validated GSE213621 merge; gene_counts.txt was the final data rename")


if __name__ == "__main__":
    main()
