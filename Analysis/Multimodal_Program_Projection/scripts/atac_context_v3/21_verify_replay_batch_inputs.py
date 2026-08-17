#!/usr/bin/env python3
"""Fail-closed preflight for one deterministic genetics replay batch."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = (
    ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=int, choices=range(1, 6), required=True)
    parser.add_argument("--candidate-root", type=Path, default=CANDIDATE)
    parser.add_argument("--fixture-mode", action="store_true")
    args = parser.parse_args()
    candidate = args.candidate_root.resolve()
    if not args.fixture_mode and candidate != CANDIDATE:
        raise RuntimeError(f"unsafe candidate root: {candidate}")
    prepared = candidate / "genetics/prepared"
    gate = read(prepared / "PROMOTED_UPSTREAM_GATE.tsv")
    if len(gate) != 1 or gate[0]["status"] != "PROMOTED":
        raise RuntimeError("promoted upstream gate is absent or invalid")
    manifest_path = prepared / "replay_input_manifest.tsv"
    if sha256(manifest_path) != gate[0]["replay_input_manifest_sha256"]:
        raise RuntimeError("replay input manifest does not match promoted gate")
    inputs = read(manifest_path)
    for row in inputs:
        path = Path(row["path"])
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256(path) != row["sha256"]
        ):
            raise RuntimeError(f"frozen replay input changed: {path}")
    batch = read(prepared / f"replay_batches/batch_{args.batch}.tsv")
    if any(int(row["batch_id"]) != args.batch for row in batch):
        raise RuntimeError("replay batch contains a different batch ID")
    keys = [(row["gwas_name"], row["ensembl"]) for row in batch]
    if len(keys) != len(set(keys)):
        raise RuntimeError("replay batch contains duplicated gene-study pairs")
    final = candidate / f"genetics/replay_execution/batch_{args.batch}"
    if final.exists() or final.is_symlink():
        raise RuntimeError(f"refusing to overwrite replay batch: {final}")
    pending = list((candidate / "genetics").glob(f".batch_{args.batch}.pending.*"))
    if pending:
        raise RuntimeError(f"unresolved prior pending replay batch: {pending}")
    print(f"REPLAY_BATCH_{args.batch}_PREFLIGHT\tPASS\tpairs={len(batch)}")


if __name__ == "__main__":
    main()
