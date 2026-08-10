#!/usr/bin/env python3
"""Extract the five validated corrected count paths for integration overrides."""

import argparse
import csv
import hashlib
import os
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def contained(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--count-manifest", required=True, type=Path)
    parser.add_argument("--baseline-overrides", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    if not (root / ".bg001_candidate_root").is_file():
        raise SystemExit("Missing candidate sentinel")
    with args.count_manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {"GSE130970", "GSE135251", "GSE174478", "GSE213621", "GSE240729"}
    if {row["dataset"] for row in rows} != expected:
        raise SystemExit("Count manifest does not contain exactly the five affected active cohorts")
    if len(rows) != len(expected):
        raise SystemExit("Count manifest contains duplicate affected cohorts")
    for row in rows:
        for path_field, hash_field in (
            ("count_path", "count_sha256"),
            ("summary_path", "summary_sha256"),
            ("log_path", "log_sha256"),
            ("environment_path", "environment_sha256"),
        ):
            path = Path(row[path_field]).resolve(strict=True)
            if not contained(path, root / "counts") or sha256(path) != row[hash_field]:
                raise SystemExit(f"Validated corrected artifact drift or path escape: {row['dataset']} {path_field}")
    with args.baseline_overrides.open(newline="") as handle:
        baseline = list(csv.DictReader(handle, delimiter="\t"))
    baseline_by_dataset = {row["dataset"]: row["count_path"] for row in baseline}
    expected_active = {
        "GSE126848", "GSE167523", "GSE135251", "GSE130970", "GSE213621",
        "GSE162694", "GSE174478", "GSE193066", "GSE240729",
    }
    if len(baseline_by_dataset) != len(baseline) or set(baseline_by_dataset) != expected_active:
        raise SystemExit("Baseline overrides do not contain exactly the nine active cohorts")
    for row in rows:
        baseline_by_dataset[row["dataset"]] = row["count_path"]
    output = args.output.resolve(strict=False)
    if not contained(output, root / "manifests") or output.parent != root / "manifests":
        raise SystemExit("Count overrides output must be directly inside RUN_ROOT/manifests")
    if output.exists() or output.is_symlink():
        raise SystemExit(f"Refusing to overwrite {args.output}")
    tmp = output.with_suffix(output.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("dataset", "count_path"), delimiter="\t")
        writer.writeheader()
        for dataset, raw_path in sorted(baseline_by_dataset.items()):
            path = Path(raw_path)
            if not path.is_file():
                raise SystemExit(f"Missing validated count matrix: {path}")
            writer.writerow({"dataset": dataset, "count_path": str(path.resolve())})
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, output)


if __name__ == "__main__":
    main()
