#!/usr/bin/env python3
"""Acquire exact author-deposited OSCAR processed inputs from a read-only commit."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.parse
import urllib.request


COMMIT = "57dce224345b3ba25baf585c8a025d4b156d5342"
FILES = (
    ("! introduction.md", 267, "c041f8dd6af4cac0fd74cc18ef2d19a9f94e0fd4"),
    ("OSCAR_DM_expression_matrix.csv.gz.00", 15728640, "7ee8a1d3fffa8953a1a04495285a73419f8c301a"),
    ("OSCAR_DM_expression_matrix.csv.gz.01", 15728640, "83cfccbda730315dafc721c883a324dadfdb31ce"),
    ("OSCAR_DM_expression_matrix.csv.gz.02", 15728640, "be9ab5ac4c862ee36781437250ca0bd14a9ee474"),
    ("OSCAR_DM_expression_matrix.csv.gz.03", 15728640, "98078b94154293277537843e9e7b7ff5207ca603"),
    ("OSCAR_DM_expression_matrix.csv.gz.04", 15728640, "fbe535f91aa752ae02b9c1716d040ad452539024"),
    ("OSCAR_DM_expression_matrix.csv.gz.05", 14142095, "e36ac68600cfd83fe88de88683ce19f9cf9d8e44"),
    ("OSCAR_EM_expression_matrix.csv.gz.00", 15728640, "0118aed7dd3f4282afcf61daf70a40d7acdcc41d"),
    ("OSCAR_EM_expression_matrix.csv.gz.01", 15728640, "f93bcbd5db9357d04fca2ded5e15e52cab739c5b"),
    ("OSCAR_EM_expression_matrix.csv.gz.02", 8538626, "9e172790908091c27cad521dddcb5fdf2141da86"),
    ("OSCAR_metadata.csv", 3651464, "cfdc8f4f9212f94a50990ae53d30efd5e7d24e26"),
)


class AcquisitionError(RuntimeError):
    """Raised when a read-only author-repository input differs."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_blob_sha1(path: Path, size: int) -> str:
    digest = hashlib.sha1(usedforsecurity=False)
    digest.update(f"blob {size}\0".encode("ascii"))
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.task_index < 0 or args.task_index >= len(FILES):
        raise AcquisitionError("OSCAR acquisition task index differs")
    if args.output.exists():
        raise AcquisitionError("OSCAR acquisition output already exists")
    name, expected_size, expected_blob = FILES[args.task_index]
    quoted_path = urllib.parse.quote(f"data/2 OSCAR/{name}")
    url = f"https://raw.githubusercontent.com/Wangxiaoyue-lab/OSCAR/{COMMIT}/{quoted_path}"
    args.output.mkdir(parents=True, exist_ok=False)
    raw_dir = args.output / "raw"
    raw_dir.mkdir()
    destination = raw_dir / name
    request = urllib.request.Request(url, headers={"User-Agent": "MASLD-model-benchmark/1"})
    with urllib.request.urlopen(request, timeout=120) as response, destination.open("wb") as handle:
        if response.status != 200 or response.geturl() != url:
            raise AcquisitionError("OSCAR immutable source response differs")
        shutil.copyfileobj(response, handle, length=8 * 1024 * 1024)
    observed_size = destination.stat().st_size
    observed_blob = git_blob_sha1(destination, observed_size)
    if observed_size != expected_size or observed_blob != expected_blob:
        raise AcquisitionError("OSCAR immutable source content differs")
    receipt = {
        "schema_version": "masld-bench-oscar-processed-input-acquisition-v1",
        "dataset_id": "cra009621_oscar",
        "repository": "Wangxiaoyue-lab/OSCAR",
        "license": "MIT",
        "commit": COMMIT,
        "repository_path": f"data/2 OSCAR/{name}",
        "url": url,
        "size_bytes": observed_size,
        "git_blob_sha1": observed_blob,
        "sha256": file_sha256(destination),
        "outcome_model_fit": False,
        "biological_replication_claim": False,
        "gse313774_accessed": False,
    }
    (args.output / "acquisition.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
