#!/usr/bin/env python3
"""Validate the five candidate matrices and freeze recount manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import stat as stat_module
import uuid
from pathlib import Path
from typing import Iterable


RUN_RE = re.compile(r"^bg001-fragment-v211-gencode49-[0-9]{8}T[0-9]{6}Z$")
EXPECTED = {
    "GSE130970": 78,
    "GSE135251": 216,
    "GSE174478": 93,
    "GSE213621": 367,
    "GSE240729": 66,
}
BAM_MANIFEST_FIELDS = (
    "dataset", "sample_id", "bam_path", "layout", "strandedness",
    "expected_status", "exclusion_reason", "bam_size", "bam_mtime", "bam_sha256",
)
CHECKSUM_FIELDS = (
    "dataset", "sample_id", "bam_path", "bam_size", "bam_mtime", "bam_sha256",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def require_confined_parent(path: Path, root: Path) -> Path:
    """Return one canonical, nonsymlinked output parent inside the candidate."""
    root = root.resolve(strict=True)
    parent = path.parent
    if not path.is_absolute() or not parent.is_dir() or parent.is_symlink():
        raise SystemExit(f"Output parent is missing, relative, or symlinked: {parent}")
    resolved = parent.resolve(strict=True)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise SystemExit(f"Output parent escapes the candidate root: {parent}") from exc
    if resolved != parent:
        raise SystemExit(f"Output parent contains symlink indirection: {parent}")
    return resolved


def require_regular_inode(path: Path, label: str) -> os.stat_result:
    try:
        observed = path.lstat()
    except FileNotFoundError as exc:
        raise SystemExit(f"{label} disappeared: {path}") from exc
    if not stat_module.S_ISREG(observed.st_mode):
        raise SystemExit(f"{label} is not a regular file: {path}")
    return observed


def same_inode(left: os.stat_result, right: os.stat_result) -> bool:
    return (left.st_dev, left.st_ino) == (right.st_dev, right.st_ino)


def atomic_write_tsv(
    path: Path,
    fieldnames: Iterable[str],
    rows: Iterable[dict[str, object]],
    *,
    root: Path,
) -> None:
    """Publish a new TSV without following or replacing any existing path."""
    parent = require_confined_parent(path, root)
    legacy_temporary = path.with_suffix(path.suffix + ".tmp")
    for forbidden in (path, legacy_temporary):
        if forbidden.exists() or forbidden.is_symlink():
            raise SystemExit(f"Refusing existing output or temporary artifact: {forbidden}")

    temporary = parent / f".{path.name}.tmp.{uuid.uuid4().hex}"
    if temporary.exists() or temporary.is_symlink():
        raise SystemExit(f"Refusing colliding UUID temporary artifact: {temporary}")
    opened: os.stat_result | None = None
    try:
        with temporary.open("x", newline="") as handle:
            opened = os.fstat(handle.fileno())
            if not stat_module.S_ISREG(opened.st_mode) or opened.st_nlink != 1:
                raise SystemExit(f"New temporary artifact is not independent: {temporary}")
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
            if not same_inode(opened, require_regular_inode(temporary, "Temporary artifact")):
                raise SystemExit(f"Temporary artifact changed while open: {temporary}")

        if require_confined_parent(path, root) != parent:
            raise SystemExit(f"Output parent changed during publication: {parent}")
        staged = require_regular_inode(temporary, "Temporary artifact")
        if opened is None or not same_inode(opened, staged) or staged.st_nlink != 1:
            raise SystemExit(f"Temporary artifact changed before publication: {temporary}")
        if path.exists() or path.is_symlink():
            raise SystemExit(f"Output appeared during publication: {path}")
        try:
            os.link(temporary, path, follow_symlinks=False)
        except FileExistsError as exc:
            raise SystemExit(f"Output appeared during publication: {path}") from exc
        fsync_directory(parent)
        published = require_regular_inode(path, "Published output")
        if not same_inode(opened, published):
            raise SystemExit(f"Published output is not the staged artifact: {path}")
        temporary.unlink()
        fsync_directory(parent)
        published = require_regular_inode(path, "Published output")
        if not same_inode(opened, published) or published.st_nlink != 1:
            raise SystemExit(f"Published output is not an independent regular file: {path}")
        if path.resolve(strict=True).parent != parent:
            raise SystemExit(f"Published output escapes its validated parent: {path}")
    except BaseException:
        # Remove only the UUID temporary that this call created.  Never unlink a
        # path that was replaced concurrently with another inode.
        if opened is not None:
            try:
                current = temporary.lstat()
            except FileNotFoundError:
                current = None
            if current is not None and same_inode(opened, current):
                temporary.unlink()
                fsync_directory(parent)
        raise


def verify_sha_manifest(path: Path, expected_names: set[str]) -> None:
    observed: set[str] = set()
    for line in path.read_text().splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or not re.fullmatch(r"[0-9a-f]{64}", fields[0]):
            raise SystemExit(f"Malformed SHA-256 manifest: {path}")
        artifact = Path(fields[1]).resolve(strict=True)
        if artifact.parent != path.parent.resolve() or artifact.name in observed:
            raise SystemExit(f"Escaping or duplicate SHA-256 entry: {artifact}")
        observed.add(artifact.name)
        if sha256(artifact) != fields[0]:
            raise SystemExit(f"Task artifact hash drift: {artifact}")
    if observed != expected_names:
        raise SystemExit(f"Task SHA-256 coverage differs in {path}: {observed}")


def environment_values(path: Path, key: str) -> list[str]:
    values = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split("\t", 1)
        if len(fields) == 2 and fields[0] == key:
            values.append(fields[1])
    if not values:
        raise SystemExit(f"Environment record lacks {key}: {path}")
    return sorted(set(values))


def included_manifest_rows(bam_rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    rows = {dataset: [] for dataset in EXPECTED}
    for row in bam_rows:
        dataset = row["dataset"]
        if dataset in rows and row["expected_status"] == "included":
            rows[dataset].append(row)
    for dataset, expected in EXPECTED.items():
        observed = rows[dataset]
        if len(observed) != expected or len({row["sample_id"] for row in observed}) != expected:
            raise SystemExit(
                f"BAM manifest sample contract differs for {dataset}: "
                f"{len(observed)} rows / {len({row['sample_id'] for row in observed})} unique"
            )
    return rows


def expected_checksum_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [{field: row[field] for field in CHECKSUM_FIELDS} for row in rows]


def task_checksum_contracts(
    root: Path,
    manifest_rows: dict[str, list[dict[str, str]]],
) -> list[tuple[Path, list[dict[str, str]]]]:
    contracts = [
        (root / "counts" / dataset, expected_checksum_rows(manifest_rows[dataset]))
        for dataset in ("GSE130970", "GSE135251", "GSE174478", "GSE240729")
    ]
    all_gse213621 = manifest_rows["GSE213621"]
    for index in range(4):
        start = index * 92
        contracts.append(
            (
                root / f"counts/GSE213621/chunks/chunk{index}",
                expected_checksum_rows(all_gse213621[start : min(start + 92, len(all_gse213621))]),
            )
        )
    return contracts


def validate_task_checksum_tables(
    root: Path,
    manifest_rows: dict[str, list[dict[str, str]]],
) -> dict[tuple[str, str], dict[str, str]]:
    """Bind every checksum row, including order, to its exact recount task."""
    observed_by_sample: dict[tuple[str, str], dict[str, str]] = {}
    for task_dir, expected_rows in task_checksum_contracts(root, manifest_rows):
        path = task_dir / "bam_checksums.tsv"
        if not path.is_file() or path.is_symlink():
            raise SystemExit(f"Missing or symlinked task BAM checksum table: {path}")
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != CHECKSUM_FIELDS:
                raise SystemExit(f"Unexpected task BAM checksum schema: {path}")
            observed_rows = list(reader)
        if observed_rows != expected_rows:
            mismatch = next(
                (
                    index
                    for index, (observed, expected) in enumerate(
                        zip(observed_rows, expected_rows), start=1
                    )
                    if observed != expected
                ),
                min(len(observed_rows), len(expected_rows)) + 1,
            )
            raise SystemExit(
                f"Task BAM checksum rows differ from the frozen manifest slice: "
                f"{path} row {mismatch}; observed={len(observed_rows)} expected={len(expected_rows)}"
            )
        for row in observed_rows:
            key = (row["dataset"], row["sample_id"])
            if key in observed_by_sample:
                raise SystemExit(f"Duplicate BAM checksum row: {key}")
            if not re.fullmatch(r"[0-9a-f]{64}", row["bam_sha256"]):
                raise SystemExit(f"Invalid BAM SHA-256: {key}")
            observed_by_sample[key] = row
    return observed_by_sample


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--artifact-origin", choices=("native", "imported"), default="native")
    parser.add_argument("--execution-source-run-id")
    parser.add_argument("--execution-source-manifest-sha256")
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    if not (root / ".bg001_candidate_root").is_file():
        raise SystemExit("Missing candidate sentinel")
    manifest_path = root / "manifests/bam_manifest.tsv"
    marker = root / "BAM_MANIFEST_FROZEN"
    if not marker.is_file() or marker.read_text().strip() != sha256(manifest_path):
        raise SystemExit("BAM manifest is not cryptographically frozen")
    contract = json.loads((root / "contract/run_contract.json").read_text())
    destination_source_manifest_sha256 = contract["source_manifest_sha256"]
    execution_source_run_id = args.execution_source_run_id or root.name
    execution_source_manifest_sha256 = (
        args.execution_source_manifest_sha256 or destination_source_manifest_sha256
    )
    if not RUN_RE.fullmatch(execution_source_run_id):
        raise SystemExit(f"Invalid execution-source run ID: {execution_source_run_id}")
    if not re.fullmatch(r"[0-9a-f]{64}", execution_source_manifest_sha256):
        raise SystemExit("Invalid execution-source manifest SHA-256")
    if args.artifact_origin == "native":
        if (
            execution_source_run_id != root.name
            or execution_source_manifest_sha256 != destination_source_manifest_sha256
        ):
            raise SystemExit("Native recount provenance must name the destination execution snapshot")
    elif (
        execution_source_run_id == root.name
        or not args.execution_source_run_id
        or not args.execution_source_manifest_sha256
    ):
        raise SystemExit("Imported recount provenance requires an explicit, different execution source")
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != BAM_MANIFEST_FIELDS:
            raise SystemExit("Unexpected BAM manifest schema")
        bam_rows = list(reader)
    manifest_rows = included_manifest_rows(bam_rows)

    expected_task_dirs = {
        root / "counts/GSE130970",
        root / "counts/GSE135251",
        root / "counts/GSE174478",
        root / "counts/GSE240729",
        *(root / f"counts/GSE213621/chunks/chunk{index}" for index in range(4)),
    }
    observed_checksum_paths = set((root / "counts").glob("**/bam_checksums.tsv"))
    expected_checksum_paths = {path / "bam_checksums.tsv" for path in expected_task_dirs}
    if observed_checksum_paths != expected_checksum_paths:
        raise SystemExit("Recount task BAM-checksum file set differs from the exact eight-task contract")
    checksum_rows = validate_task_checksum_tables(root, manifest_rows)
    task_environment_hashes: set[str] = set()
    task_featurecounts_hashes: set[str] = set()
    task_samtools_hashes: set[str] = set()
    task_samtools_versions: set[str] = set()
    for task_dir in sorted(expected_task_dirs):
        complete = task_dir / "COMPLETE"
        provenance = task_dir / "provenance.sha256"
        if not complete.is_file() or complete.is_symlink() or not provenance.is_file():
            raise SystemExit(f"Task lacks COMPLETE/provenance: {task_dir}")
        verify_sha_manifest(
            provenance,
            {"bam_checksums.tsv", "environment.txt", "gene_counts.txt", "gene_counts.txt.summary", "featureCounts.log", "validation.json"},
        )
        environment = task_dir / "environment.txt"
        task_environment_hashes.add(sha256(environment))
        task_featurecounts_hashes.update(environment_values(environment, "featureCounts_executable_sha256"))
        task_samtools_hashes.update(environment_values(environment, "samtools_executable_sha256"))
        task_samtools_versions.update(environment_values(environment, "samtools_version"))
    expected_checksum_count = sum(EXPECTED.values())
    if len(checksum_rows) != expected_checksum_count:
        raise SystemExit(
            f"Expected {expected_checksum_count} BAM checksums, found {len(checksum_rows)}"
        )
    if not all(
        len(values) == 1
        for values in (
            task_environment_hashes,
            task_featurecounts_hashes,
            task_samtools_hashes,
            task_samtools_versions,
        )
    ):
        raise SystemExit("The eight recount tasks did not use one identical pinned tool environment")
    pinned_environment_sha256 = next(iter(task_environment_hashes))
    pinned_featurecounts_sha256 = next(iter(task_featurecounts_hashes))
    pinned_samtools_sha256 = next(iter(task_samtools_hashes))
    pinned_samtools_version = next(iter(task_samtools_versions))

    for row in bam_rows:
        if row["expected_status"] != "included":
            continue
        key = (row["dataset"], row["sample_id"])
        check = checksum_rows.get(key)
        if not check:
            raise SystemExit(f"Missing BAM checksum: {key}")
        if (
            check["bam_path"] != row["bam_path"]
            or check["bam_size"] != row["bam_size"]
            or check["bam_mtime"] != row["bam_mtime"]
            or check["bam_sha256"] != row["bam_sha256"]
        ):
            raise SystemExit(f"BAM provenance mismatch: {key}")

    count_rows = []
    for dataset, expected_samples in EXPECTED.items():
        dataset_dir = root / "counts" / dataset
        validation_path = dataset_dir / "validation.json"
        completion = dataset_dir / ("MERGE_COMPLETE" if dataset == "GSE213621" else "COMPLETE")
        if not validation_path.is_file() or not completion.is_file():
            raise SystemExit(f"Missing final validation for {dataset}: {validation_path}")
        validation = json.loads(validation_path.read_text())
        if validation.get("status") != "PASS" or validation.get("n_samples") != expected_samples:
            raise SystemExit(f"Invalid final validation for {dataset}")
        counts = dataset_dir / "gene_counts.txt"
        summary = dataset_dir / "gene_counts.txt.summary"
        log = dataset_dir / "featureCounts.log"
        environment = dataset_dir / "environment.txt"
        if (
            validation.get("counts_sha256") != sha256(counts)
            or validation.get("summary_sha256") != sha256(summary)
            or validation.get("log_sha256") != sha256(log)
            or validation.get("environment_sha256") != sha256(environment)
            or validation.get("n_genes") != 86_369
            or validation.get("featurecounts_version") != "2.1.1"
            or validation.get("gtf_sha256") != "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
        ):
            raise SystemExit(f"Validated output hash drift for {dataset}")
        if dataset == "GSE213621":
            source_commands = []
            for chunk_index in range(4):
                chunk_validation = dataset_dir / "chunks" / f"chunk{chunk_index}" / "validation.json"
                payload = json.loads(chunk_validation.read_text())
                source_commands.append(payload["command"])
        else:
            source_commands = [validation["command"]]
        count_rows.append(
            {
                "dataset": dataset,
                "count_path": str(counts),
                "summary_path": str(summary),
                "log_path": str(log),
                "environment_path": str(environment),
                "n_samples": expected_samples,
                "n_genes": validation["n_genes"],
                "featurecounts_version": "2.1.1",
                "command_contract": "-p --countReadPairs -B -s 2",
                "exact_command_json": json.dumps(source_commands, separators=(",", ":")),
                "gtf_sha256": "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
                "bam_manifest_sha256": sha256(manifest_path),
                # This legacy-named field remains the manifest of the snapshot
                # that actually executed featureCounts.  It must never be
                # rewritten to imply that imported counts were produced here.
                "source_manifest_sha256": execution_source_manifest_sha256,
                "artifact_origin": args.artifact_origin,
                "execution_source_run_id": execution_source_run_id,
                "execution_source_manifest_sha256": execution_source_manifest_sha256,
                "destination_analysis_source_manifest_sha256": destination_source_manifest_sha256,
                "pinned_environment_sha256": pinned_environment_sha256,
                "pinned_featurecounts_executable_sha256": pinned_featurecounts_sha256,
                "pinned_samtools_executable_sha256": pinned_samtools_sha256,
                "pinned_samtools_version": pinned_samtools_version,
                "featurecounts_executable_sha256_json": json.dumps(environment_values(environment, "featureCounts_executable_sha256"), separators=(",", ":")),
                "samtools_version_json": json.dumps(environment_values(environment, "samtools_version"), separators=(",", ":")),
                "samtools_executable_sha256_json": json.dumps(environment_values(environment, "samtools_executable_sha256"), separators=(",", ":")),
                "count_sha256": sha256(counts),
                "summary_sha256": sha256(summary),
                "log_sha256": sha256(log),
                "environment_sha256": sha256(environment),
                "validation_sha256": sha256(validation_path),
            }
        )
    count_manifest = root / "manifests/count_manifest.tsv"
    atomic_write_tsv(
        count_manifest,
        count_rows[0].keys(),
        count_rows,
        root=root,
    )
    print(f"PASS: validated eight COMPLETE recount tasks and froze {count_manifest}")


if __name__ == "__main__":
    main()
