#!/usr/bin/env python3
"""Verify candidate containment and every frozen source byte before execution."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


RUN_RE = re.compile(r"^bg001-fragment-v211-gencode49-[0-9]{8}T[0-9]{6}Z$")
MANIFEST_FIELDS = (
    "dataset", "sample_id", "bam_path", "layout", "strandedness",
    "expected_status", "exclusion_reason", "bam_size", "bam_mtime", "bam_sha256",
)
EXPECTED = {
    "GSE130970": 78,
    "GSE135251": 216,
    "GSE174478": 93,
    "GSE213621": 367,
    "GSE240729": 66,
}
DOCUMENTED_MISSING = {
    ("GSE174478", "SRR14551000"),
    ("GSE240729", "SRR25630203"),
}


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


def validate_bam_manifest(path: Path, require_hashes: bool) -> list[dict[str, str]]:
    if b"\r" in path.read_bytes():
        raise SystemExit(f"BAM manifest must use LF-only line endings: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise SystemExit("BAM manifest schema/order differs from the contract")
        rows = list(reader)
    if len(rows) != sum(EXPECTED.values()) + len(DOCUMENTED_MISSING):
        raise SystemExit("BAM manifest row count differs from the 820+2 contract")
    included = [row for row in rows if row["expected_status"] == "included"]
    unavailable = [row for row in rows if row["expected_status"] == "documented_unavailable"]
    if {(row["dataset"], row["sample_id"]) for row in unavailable} != DOCUMENTED_MISSING:
        raise SystemExit("Documented unavailable BAM rows differ from the contract")
    if any(row["bam_path"] or row["bam_size"] or row["bam_mtime"] or row["bam_sha256"] for row in unavailable):
        raise SystemExit("Documented unavailable rows must not contain BAM provenance")
    keys = [(row["dataset"], row["sample_id"]) for row in included]
    paths = [row["bam_path"] for row in included]
    if len(keys) != len(set(keys)) or len(paths) != len(set(paths)):
        raise SystemExit("Included BAM manifest rows contain duplicate samples or paths")
    for dataset, expected in EXPECTED.items():
        if sum(row["dataset"] == dataset for row in included) != expected:
            raise SystemExit(f"{dataset} BAM manifest count differs from {expected}")
    for row in included:
        if row["layout"] != "paired" or row["strandedness"] != "2" or row["exclusion_reason"]:
            raise SystemExit(f"Invalid included BAM semantics: {row['dataset']}/{row['sample_id']}")
        bam = Path(row["bam_path"])
        if not bam.is_absolute() or not bam.is_file() or bam.is_symlink():
            raise SystemExit(f"Included BAM is missing, relative, or symlinked: {bam}")
        stat = bam.stat()
        if row["bam_size"] != str(stat.st_size) or row["bam_mtime"] != str(int(stat.st_mtime)):
            raise SystemExit(f"Included BAM size/mtime drift: {bam}")
        digest = row["bam_sha256"]
        if require_hashes and not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise SystemExit(f"Included BAM lacks frozen SHA-256: {bam}")
        if not require_hashes and digest:
            raise SystemExit(f"Draft BAM manifest unexpectedly contains SHA-256: {bam}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--require-bam-frozen", action="store_true")
    args = parser.parse_args()

    root = args.run_root.resolve(strict=True)
    if (root / "REJECTED").exists():
        raise SystemExit("Candidate is permanently marked REJECTED")
    sentinel = root / ".bg001_candidate_root"
    if not sentinel.is_file() or sentinel.is_symlink():
        raise SystemExit("Candidate sentinel is missing or is a symlink")
    run_id = sentinel.read_text().strip()
    if not RUN_RE.fullmatch(run_id) or root.name != run_id:
        raise SystemExit("Candidate run ID and sentinel disagree")

    contract_path = root / "contract/run_contract.json"
    contract = json.loads(contract_path.read_text())
    if Path(contract["candidate_root"]).resolve(strict=True) != root:
        raise SystemExit("Contract candidate_root does not resolve to this run")
    project = Path(contract["project_root"]).resolve(strict=True)
    allowed_parent = (project / "results/remediation/bg001").resolve(strict=True)
    if not contained(root, allowed_parent) or root == allowed_parent:
        raise SystemExit("Candidate root escapes the fixed BG-001 remediation parent")
    if sha256(root / "contract/dirty_worktree.patch") != contract["dirty_patch_sha256"]:
        raise SystemExit("Authoritative dirty-worktree patch hash drifted")
    if sha256(root / "contract/source_manifest.tsv") != contract["source_manifest_sha256"]:
        raise SystemExit("Frozen source-manifest hash drifted")
    if sha256(root / "contract/source_regression_check.txt") != contract["source_regression_check_sha256"]:
        raise SystemExit("Frozen producer regression-check record drifted")
    if sha256(root / "contract/dependency_source_inventory.tsv") != contract["dependency_source_inventory_sha256"]:
        raise SystemExit("Frozen dependency inventory hash drifted")
    if sha256(root / "contract/untracked_dependency_sources.tsv") != contract["untracked_dependency_sources_sha256"]:
        raise SystemExit("Frozen untracked-source inventory hash drifted")

    draft = root / "manifests/bam_manifest.draft.tsv"
    if not draft.is_file() or sha256(draft) != contract["bam_manifest_draft_sha256"]:
        raise SystemExit("Draft BAM manifest is missing or changed")
    validate_bam_manifest(draft, require_hashes=False)
    frozen_marker = root / "BAM_MANIFEST_FROZEN"
    final_manifest = root / "manifests/bam_manifest.tsv"
    freeze_record = root / "manifests/bam_manifest.freeze.json"
    if frozen_marker.exists():
        if frozen_marker.is_symlink() or not final_manifest.is_file() or not freeze_record.is_file():
            raise SystemExit("Frozen BAM-manifest transaction is incomplete")
        final_sha = sha256(final_manifest)
        if frozen_marker.read_text().strip() != final_sha:
            raise SystemExit("Frozen BAM-manifest marker digest differs")
        record = json.loads(freeze_record.read_text())
        if (
            record.get("draft_sha256") != contract["bam_manifest_draft_sha256"]
            or record.get("manifest_sha256") != final_sha
            or record.get("included_bams") != sum(EXPECTED.values())
        ):
            raise SystemExit("Frozen BAM-manifest record differs from the run contract")
        shard_hashes = record.get("shard_sha256", {})
        expected_shards = {f"task_{index}.tsv" for index in range(8)}
        if set(shard_hashes) != expected_shards:
            raise SystemExit("Frozen BAM-manifest record does not name exactly eight hash shards")
        for name, expected_hash in shard_hashes.items():
            shard = root / "manifests/bam_hash_shards" / name
            if not shard.is_file() or sha256(shard) != expected_hash:
                raise SystemExit(f"Frozen BAM hash shard drift: {name}")
        validate_bam_manifest(final_manifest, require_hashes=True)
    elif final_manifest.exists() or freeze_record.exists() or args.require_bam_frozen:
        raise SystemExit("Final BAM manifest is absent or was published without its commit marker")

    # Candidate artifacts must never redirect writes through a symlink.
    links = [path for path in root.rglob("*") if path.is_symlink()]
    if links:
        raise SystemExit(f"Candidate tree contains symlinks: {links[:5]}")

    snapshot = root / "source_snapshot"
    manifest_path = root / "contract/source_manifest.tsv"
    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {row["relative_path"] for row in rows}
    observed = {
        str(path.relative_to(snapshot))
        for path in snapshot.rglob("*")
        if path.is_file()
    }
    if observed != expected:
        raise SystemExit(
            f"Frozen source file set drifted: missing={sorted(expected-observed)[:5]}, "
            f"extra={sorted(observed-expected)[:5]}"
        )
    for row in rows:
        path = snapshot / row["relative_path"]
        if path.stat().st_size != int(row["size_bytes"]):
            raise SystemExit(f"Frozen source size drift: {row['relative_path']}")
        if sha256(path) != row["sha256"]:
            raise SystemExit(f"Frozen source hash drift: {row['relative_path']}")
    bam_state = "frozen" if frozen_marker.exists() else "draft"
    print(f"PASS: {run_id}; {len(rows)} frozen source files verified; BAM manifest {bam_state}")


if __name__ == "__main__":
    main()
