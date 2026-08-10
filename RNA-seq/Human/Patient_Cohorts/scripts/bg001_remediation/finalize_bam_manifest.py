#!/usr/bin/env python3
"""Atomically publish the cryptographically frozen BAM execution manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path


EXPECTED = {
    "GSE130970": 78,
    "GSE135251": 216,
    "GSE174478": 93,
    "GSE213621": 367,
    "GSE240729": 66,
}
SHARD_FIELDS = ("task", "dataset", "sample_id", "bam_path", "bam_size", "bam_mtime", "bam_sha256")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_task_rows(rows: list[dict[str, str]], task: int) -> list[dict[str, str]]:
    if task <= 3:
        dataset = ("GSE130970", "GSE135251", "GSE174478", "GSE240729")[task]
        return [row for row in rows if row["dataset"] == dataset and row["expected_status"] == "included"]
    dataset_rows = [row for row in rows if row["dataset"] == "GSE213621" and row["expected_status"] == "included"]
    start = (task - 4) * 92
    return dataset_rows[start : min(start + 92, len(dataset_rows))]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    if not (root / ".bg001_candidate_root").is_file():
        raise SystemExit("Missing candidate sentinel")
    draft = root / "manifests/bam_manifest.draft.tsv"
    contract = json.loads((root / "contract/run_contract.json").read_text())
    if sha256(draft) != contract["bam_manifest_draft_sha256"]:
        raise SystemExit("Draft BAM manifest hash differs from the run contract")
    with draft.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        rows = list(reader)

    shard_dir = root / "manifests/bam_hash_shards"
    expected_paths = {shard_dir / f"task_{task}.tsv" for task in range(8)}
    observed_paths = set(shard_dir.glob("*"))
    if observed_paths != expected_paths:
        raise SystemExit(
            f"BAM hash shards differ: missing={sorted(str(x) for x in expected_paths-observed_paths)}, "
            f"extra={sorted(str(x) for x in observed_paths-expected_paths)}"
        )

    hash_by_key: dict[tuple[str, str], str] = {}
    shard_hashes: dict[str, str] = {}
    for task in range(8):
        path = shard_dir / f"task_{task}.tsv"
        shard_hashes[path.name] = sha256(path)
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != SHARD_FIELDS:
                raise SystemExit(f"BAM hash shard schema differs: {path}")
            observed = list(reader)
        expected_rows = expected_task_rows(rows, task)
        if len(observed) != len(expected_rows):
            raise SystemExit(f"Task {task} hash cardinality differs")
        for expected_row, observed_row in zip(expected_rows, observed):
            key = (expected_row["dataset"], expected_row["sample_id"])
            if key in hash_by_key:
                raise SystemExit(f"Duplicate BAM hash row: {key}")
            for field in ("dataset", "sample_id", "bam_path", "bam_size", "bam_mtime"):
                if observed_row[field] != expected_row[field]:
                    raise SystemExit(f"Task {task} BAM provenance differs for {key}: {field}")
            if observed_row["task"] != str(task) or not re.fullmatch(r"[0-9a-f]{64}", observed_row["bam_sha256"]):
                raise SystemExit(f"Task {task} has invalid task or SHA-256 for {key}")
            bam = Path(expected_row["bam_path"])
            stat = bam.stat()
            if str(stat.st_size) != expected_row["bam_size"] or str(int(stat.st_mtime)) != expected_row["bam_mtime"]:
                raise SystemExit(f"BAM changed after prehash: {bam}")
            hash_by_key[key] = observed_row["bam_sha256"]

    if len(hash_by_key) != sum(EXPECTED.values()):
        raise SystemExit(f"Expected 820 unique BAM hashes, found {len(hash_by_key)}")
    for row in rows:
        if row["expected_status"] == "included":
            row["bam_sha256"] = hash_by_key[(row["dataset"], row["sample_id"])]

    manifest = root / "manifests/bam_manifest.tsv"
    record_path = root / "manifests/bam_manifest.freeze.json"
    marker = root / "BAM_MANIFEST_FROZEN"
    for path in (manifest, record_path, marker):
        if path.exists() or path.is_symlink():
            raise SystemExit(f"Refusing existing BAM-manifest publication artifact: {path}")
    manifest_tmp = manifest.with_suffix(".tsv.tmp")
    with manifest_tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    manifest_sha = sha256(manifest_tmp)
    record = {
        "schema_version": "1.0",
        "run_id": root.name,
        "draft_sha256": sha256(draft),
        "manifest_sha256": manifest_sha,
        "included_bams": len(hash_by_key),
        "documented_unavailable": len(rows) - len(hash_by_key),
        "shard_sha256": shard_hashes,
        "finalized_utc": datetime.now(timezone.utc).isoformat(),
    }
    record_tmp = record_path.with_suffix(".json.tmp")
    record_tmp.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    os.replace(manifest_tmp, manifest)
    os.replace(record_tmp, record_path)
    marker_tmp = marker.with_suffix(".tmp")
    marker_tmp.write_text(manifest_sha + "\n")
    os.replace(marker_tmp, marker)
    print(f"PASS: atomically froze 820 BAM SHA-256 values in {manifest}")


if __name__ == "__main__":
    main()
