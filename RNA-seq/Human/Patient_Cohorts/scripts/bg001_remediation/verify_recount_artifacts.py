#!/usr/bin/env python3
"""Seal or reverify every recount artifact consumed by the four analysis arms."""

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


AFFECTED = {"GSE130970", "GSE135251", "GSE174478", "GSE213621", "GSE240729"}
EXPECTED_SAMPLES = {
    "GSE130970": 78,
    "GSE135251": 216,
    "GSE174478": 93,
    "GSE213621": 367,
    "GSE240729": 66,
}
EXPECTED_GTF_SHA256 = "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
RUN_RE = re.compile(r"^bg001-fragment-v211-gencode49-[0-9]{8}T[0-9]{6}Z$")
LEGACY_COUNT_MANIFEST_FIELDS = (
    "dataset", "count_path", "summary_path", "log_path", "environment_path",
    "n_samples", "n_genes", "featurecounts_version", "command_contract",
    "exact_command_json", "gtf_sha256", "bam_manifest_sha256",
    "source_manifest_sha256", "pinned_environment_sha256",
    "pinned_featurecounts_executable_sha256", "pinned_samtools_executable_sha256",
    "pinned_samtools_version", "featurecounts_executable_sha256_json",
    "samtools_version_json", "samtools_executable_sha256_json", "count_sha256",
    "summary_sha256", "log_sha256", "environment_sha256", "validation_sha256",
)
COUNT_MANIFEST_FIELDS = (
    "dataset", "count_path", "summary_path", "log_path", "environment_path",
    "n_samples", "n_genes", "featurecounts_version", "command_contract",
    "exact_command_json", "gtf_sha256", "bam_manifest_sha256",
    "source_manifest_sha256", "artifact_origin", "execution_source_run_id",
    "execution_source_manifest_sha256", "destination_analysis_source_manifest_sha256",
    "pinned_environment_sha256", "pinned_featurecounts_executable_sha256",
    "pinned_samtools_executable_sha256", "pinned_samtools_version",
    "featurecounts_executable_sha256_json", "samtools_version_json",
    "samtools_executable_sha256_json", "count_sha256", "summary_sha256",
    "log_sha256", "environment_sha256", "validation_sha256",
)
ACTIVE = {
    "GSE126848", "GSE167523", "GSE135251", "GSE130970", "GSE213621",
    "GSE162694", "GSE174478", "GSE193066", "GSE240729",
}
DIRECT = {"GSE130970", "GSE135251", "GSE174478", "GSE240729"}
TASK_FILES = {
    "bam_checksums.tsv", "environment.txt", "gene_counts.txt", "gene_counts.txt.summary",
    "featureCounts.log", "validation.json", "provenance.sha256", "COMPLETE",
}
MERGED_FILES = {
    "gene_counts.txt", "gene_counts.txt.summary", "featureCounts.log", "environment.txt",
    "validation.json", "merge_publication.json", "MERGE_COMPLETE",
}
BAM_MANIFEST_FIELDS = (
    "dataset", "sample_id", "bam_path", "layout", "strandedness",
    "expected_status", "exclusion_reason", "bam_size", "bam_mtime", "bam_sha256",
)
CHECKSUM_FIELDS = (
    "dataset", "sample_id", "bam_path", "bam_size", "bam_mtime", "bam_sha256",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def require_new_output_target(path: Path, root: Path, legacy_temporary: Path) -> Path:
    parent = require_confined_parent(path, root)
    for forbidden in (path, legacy_temporary):
        if forbidden.exists() or forbidden.is_symlink():
            raise SystemExit(f"Refusing existing output or temporary artifact: {forbidden}")
    return parent


def atomic_write_text(
    path: Path,
    text: str,
    *,
    root: Path,
    legacy_temporary: Path,
) -> None:
    """Publish a new text artifact without following or replacing existing paths."""
    parent = require_new_output_target(path, root, legacy_temporary)
    temporary = parent / f".{path.name}.tmp.{uuid.uuid4().hex}"
    if temporary.exists() or temporary.is_symlink():
        raise SystemExit(f"Refusing colliding UUID temporary artifact: {temporary}")
    opened: os.stat_result | None = None
    try:
        with temporary.open("x") as handle:
            opened = os.fstat(handle.fileno())
            if not stat_module.S_ISREG(opened.st_mode) or opened.st_nlink != 1:
                raise SystemExit(f"New temporary artifact is not independent: {temporary}")
            handle.write(text)
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
        if path.resolve(strict=True).parent != parent or sha256(path) != hashlib.sha256(text.encode()).hexdigest():
            raise SystemExit(f"Published output differs or escapes its validated parent: {path}")
    except BaseException:
        if opened is not None:
            try:
                current = temporary.lstat()
            except FileNotFoundError:
                current = None
            if current is not None and same_inode(opened, current):
                temporary.unlink()
                fsync_directory(parent)
        raise


def environment_values(path: Path, key: str) -> set[str]:
    values: set[str] = set()
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split("\t", 1)
        if len(fields) == 2 and fields[0] == key:
            values.add(fields[1])
    return values


def relative_files(root: Path) -> set[str]:
    return {str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()}


def expected_count_files(root: Path) -> set[Path]:
    paths = {root / "counts" / dataset / name for dataset in DIRECT for name in TASK_FILES}
    for index in range(4):
        paths.update(root / f"counts/GSE213621/chunks/chunk{index}" / name for name in TASK_FILES)
    paths.update(root / "counts/GSE213621" / name for name in MERGED_FILES)
    return paths


def authoritative_manifest_files(root: Path) -> set[Path]:
    paths = {
        root / "manifests/bam_manifest.draft.tsv",
        root / "manifests/bam_manifest.tsv",
        root / "manifests/bam_manifest.freeze.json",
        root / "manifests/bam_manifest_finalize.preflight.json",
        root / "manifests/bam_manifest_finalize.postflight.json",
        root / "manifests/count_manifest.tsv",
        root / "manifests/count_overrides.tsv",
        root / "BAM_MANIFEST_FROZEN",
    }
    paths.update(root / f"manifests/bam_hash_shards/task_{index}.tsv" for index in range(8))
    return paths


def reject_recount_symlinks(root: Path) -> None:
    """Reject indirection anywhere in the two authoritative recount trees."""
    for relative in ("counts", "manifests"):
        base = root / relative
        if not base.is_dir() or base.is_symlink():
            raise SystemExit(f"Missing or symlinked recount directory: {base}")
        symlinks = sorted(path for path in base.rglob("*") if path.is_symlink())
        if symlinks:
            raise SystemExit(f"Recount tree contains symlink: {symlinks[0]}")


def included_manifest_rows(path: Path) -> dict[str, list[dict[str, str]]]:
    rows: dict[str, list[dict[str, str]]] = {dataset: [] for dataset in AFFECTED}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != BAM_MANIFEST_FIELDS:
            raise SystemExit("Unexpected BAM manifest schema")
        for row in reader:
            dataset = row["dataset"]
            if dataset in rows and row["expected_status"] == "included":
                rows[dataset].append(row)
    for dataset, expected_count in EXPECTED_SAMPLES.items():
        observed = rows[dataset]
        unique_samples = {row["sample_id"] for row in observed}
        if len(observed) != expected_count or len(unique_samples) != expected_count:
            raise SystemExit(
                f"BAM manifest sample contract differs for {dataset}: "
                f"{len(observed)} rows / {len(unique_samples)} unique"
            )
    return rows


def included_manifest_samples(path: Path) -> dict[str, list[str]]:
    return {
        dataset: [row["sample_id"] for row in rows]
        for dataset, rows in included_manifest_rows(path).items()
    }


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
) -> None:
    """Re-derive the exact task-level checksum attribution from the BAM manifest."""
    observed_keys: set[tuple[str, str]] = set()
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
            if key in observed_keys or not re.fullmatch(r"[0-9a-f]{64}", row["bam_sha256"]):
                raise SystemExit(f"Duplicate or invalid task BAM checksum row: {key}")
            observed_keys.add(key)
    if len(observed_keys) != sum(EXPECTED_SAMPLES.values()):
        raise SystemExit(f"Task BAM checksum coverage differs: {len(observed_keys)}")


def read_command_header(path: Path) -> str:
    with path.open(errors="replace") as handle:
        return handle.readline().rstrip("\n")


def validate_validation_payload(
    *,
    validation_path: Path,
    dataset: str,
    expected_samples: list[str],
    counts: Path,
    summary: Path,
    log: Path,
    environment: Path,
    expected_command: str,
) -> dict[str, object]:
    try:
        payload = json.loads(validation_path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        raise SystemExit(f"Unreadable validation payload: {validation_path}") from exc
    if not isinstance(payload, dict):
        raise SystemExit(f"Validation payload is not an object: {validation_path}")
    required = {
        "status": "PASS",
        "dataset": dataset,
        "n_samples": len(expected_samples),
        "n_genes": 86_369,
        "featurecounts_version": "2.1.1",
        "gtf_sha256": EXPECTED_GTF_SHA256,
        "samples": expected_samples,
        "command": expected_command,
        "counts_sha256": sha256(counts),
        "summary_sha256": sha256(summary),
        "log_sha256": sha256(log),
        "environment_sha256": sha256(environment),
    }
    for field, expected in required.items():
        if payload.get(field) != expected:
            raise SystemExit(
                f"Validation payload differs for {dataset}/{validation_path.parent.name}/"
                f"{validation_path.name}: {field}"
            )
    if read_command_header(counts) != expected_command:
        raise SystemExit(f"Count header and validation command differ: {counts}")
    return payload


def validate_dataset_artifacts(
    root: Path,
    dataset: str,
    row: dict[str, str],
    commands: list[str],
    expected_samples: list[str],
) -> None:
    """Bind one manifest row to its exact cohort files and validation evidence."""
    dataset_root = root / "counts" / dataset
    expected_paths = {
        "count_path": dataset_root / "gene_counts.txt",
        "summary_path": dataset_root / "gene_counts.txt.summary",
        "log_path": dataset_root / "featureCounts.log",
        "environment_path": dataset_root / "environment.txt",
    }
    digest_fields = {
        "count_path": "count_sha256",
        "summary_path": "summary_sha256",
        "log_path": "log_sha256",
        "environment_path": "environment_sha256",
    }
    for field, expected in expected_paths.items():
        raw_path = Path(row[field])
        if not raw_path.is_absolute() or raw_path != expected:
            raise SystemExit(
                f"Count manifest path is not the exact {dataset} destination: {field}={raw_path}"
            )
        if not raw_path.is_file() or raw_path.is_symlink():
            raise SystemExit(f"Missing or symlinked count manifest artifact: {raw_path}")
        if sha256(raw_path) != row[digest_fields[field]]:
            raise SystemExit(f"Count manifest hash drift: {dataset}/{field}")

    validation_path = dataset_root / "validation.json"
    if not validation_path.is_file() or validation_path.is_symlink():
        raise SystemExit(f"Missing or symlinked validation payload: {validation_path}")
    if sha256(validation_path) != row["validation_sha256"]:
        raise SystemExit(f"Validation hash drift: {dataset}")

    if dataset != "GSE213621":
        validate_validation_payload(
            validation_path=validation_path,
            dataset=dataset,
            expected_samples=expected_samples,
            counts=expected_paths["count_path"],
            summary=expected_paths["summary_path"],
            log=expected_paths["log_path"],
            environment=expected_paths["environment_path"],
            expected_command=commands[0],
        )
        return

    merged_command = read_command_header(expected_paths["count_path"])
    if (
        not merged_command.startswith("# Program:featureCounts v2.1.1; BG-001 keyed merge;")
        or "chunks=4" not in merged_command
    ):
        raise SystemExit("GSE213621 merged command does not identify the keyed four-chunk merge")
    for token in ('"-p"', '"--countReadPairs"', '"-B"', '"-s" "2"'):
        if token not in merged_command:
            raise SystemExit(f"GSE213621 merged command lacks {token}")
    validate_validation_payload(
        validation_path=validation_path,
        dataset=dataset,
        expected_samples=expected_samples,
        counts=expected_paths["count_path"],
        summary=expected_paths["summary_path"],
        log=expected_paths["log_path"],
        environment=expected_paths["environment_path"],
        expected_command=merged_command,
    )
    for chunk_index, command in enumerate(commands):
        start = chunk_index * 92
        chunk_samples = expected_samples[start : min(start + 92, len(expected_samples))]
        chunk_root = dataset_root / "chunks" / f"chunk{chunk_index}"
        validate_validation_payload(
            validation_path=chunk_root / "validation.json",
            dataset=dataset,
            expected_samples=chunk_samples,
            counts=chunk_root / "gene_counts.txt",
            summary=chunk_root / "gene_counts.txt.summary",
            log=chunk_root / "featureCounts.log",
            environment=chunk_root / "environment.txt",
            expected_command=command,
        )


def validate_semantics(
    root: Path,
    *,
    allow_pending_import: bool = False,
) -> tuple[set[Path], dict[str, object]]:
    root = root.resolve(strict=True)
    if not (root / ".bg001_candidate_root").is_file():
        raise SystemExit("Missing candidate sentinel")
    from bam_manifest_finalization_guard import GuardError, verify_artifacts

    try:
        verify_artifacts(root)
    except GuardError as exc:
        raise SystemExit(f"BAM-manifest guard transaction is invalid: {exc}") from exc
    reject_recount_symlinks(root)
    expected_counts = expected_count_files(root)
    observed_counts = {path for path in (root / "counts").rglob("*") if path.is_file()}
    if observed_counts != expected_counts:
        raise SystemExit(
            f"Count artifact set differs: missing={sorted(str(x) for x in expected_counts-observed_counts)[:5]}, "
            f"extra={sorted(str(x) for x in observed_counts-expected_counts)[:5]}"
        )
    manifests = authoritative_manifest_files(root)
    missing_manifests = [path for path in manifests if not path.is_file() or path.is_symlink()]
    if missing_manifests:
        raise SystemExit(f"Missing/symlinked authoritative recount manifests: {missing_manifests[:5]}")
    if any(path.name == "FAILED" or ".tmp" in path.name for path in observed_counts):
        raise SystemExit("Recount tree contains failure or temporary artifacts")

    bam_manifest = root / "manifests/bam_manifest.tsv"
    if (root / "BAM_MANIFEST_FROZEN").read_text().strip() != sha256(bam_manifest):
        raise SystemExit("BAM manifest commit marker differs")
    manifest_rows = included_manifest_rows(bam_manifest)
    manifest_samples = {
        dataset: [row["sample_id"] for row in rows]
        for dataset, rows in manifest_rows.items()
    }
    validate_task_checksum_tables(root, manifest_rows)
    with (root / "manifests/count_manifest.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        manifest_fields = tuple(reader.fieldnames or ())
        if manifest_fields not in (LEGACY_COUNT_MANIFEST_FIELDS, COUNT_MANIFEST_FIELDS):
            raise SystemExit("Count manifest schema/order differs from the frozen contract")
        count_rows = list(reader)
    if len(count_rows) != 5 or {row["dataset"] for row in count_rows} != AFFECTED:
        raise SystemExit("Count manifest does not contain exactly five affected cohorts")
    by_dataset = {row["dataset"]: row for row in count_rows}
    contract = json.loads((root / "contract/run_contract.json").read_text())
    if manifest_fields == LEGACY_COUNT_MANIFEST_FIELDS:
        artifact_origins = {"native"}
        execution_source_run_ids = {root.name}
        execution_source_manifest_hashes = {row["source_manifest_sha256"] for row in count_rows}
        destination_analysis_manifest_hashes = set(execution_source_manifest_hashes)
    else:
        artifact_origins = {row["artifact_origin"] for row in count_rows}
        execution_source_run_ids = {row["execution_source_run_id"] for row in count_rows}
        execution_source_manifest_hashes = {
            row["execution_source_manifest_sha256"] for row in count_rows
        }
        destination_analysis_manifest_hashes = {
            row["destination_analysis_source_manifest_sha256"] for row in count_rows
        }
        if any(
            row["source_manifest_sha256"] != row["execution_source_manifest_sha256"]
            for row in count_rows
        ):
            raise SystemExit("Count manifest legacy source digest is not the execution-source digest")
    if (
        len(artifact_origins) != 1
        or artifact_origins - {"native", "imported"}
        or len(execution_source_run_ids) != 1
        or not RUN_RE.fullmatch(next(iter(execution_source_run_ids), ""))
        or len(execution_source_manifest_hashes) != 1
        or len(destination_analysis_manifest_hashes) != 1
        or any(
            not re.fullmatch(r"[0-9a-f]{64}", value)
            for value in execution_source_manifest_hashes | destination_analysis_manifest_hashes
        )
    ):
        raise SystemExit("Count manifest does not name one coherent execution/analysis provenance")
    artifact_origin = next(iter(artifact_origins))
    execution_source_run_id = next(iter(execution_source_run_ids))
    execution_source_manifest_sha256 = next(iter(execution_source_manifest_hashes))
    destination_analysis_source_manifest_sha256 = next(iter(destination_analysis_manifest_hashes))
    if destination_analysis_source_manifest_sha256 != contract["source_manifest_sha256"]:
        raise SystemExit("Count manifest destination analysis snapshot differs from this run contract")
    if artifact_origin == "native":
        if (
            execution_source_run_id != root.name
            or execution_source_manifest_sha256 != contract["source_manifest_sha256"]
        ):
            raise SystemExit("Native count manifest does not name this run as its execution source")
    else:
        if execution_source_run_id == root.name:
            raise SystemExit("Imported count manifest cannot name the destination as its execution source")
        imported_marker = root / "RECOUNT_IMPORTED"
        if allow_pending_import:
            if not (root / "RECOUNT_IMPORT_IN_PROGRESS").is_file():
                raise SystemExit("Pending imported validation lacks RECOUNT_IMPORT_IN_PROGRESS")
        elif not imported_marker.is_file() or imported_marker.is_symlink():
            raise SystemExit("Imported count manifest lacks a regular RECOUNT_IMPORTED marker")
        if imported_marker.is_file() and not imported_marker.is_symlink():
            imported = json.loads(imported_marker.read_text())
            if (
                imported.get("execution_source_run_id") != execution_source_run_id
                or imported.get("execution_source_manifest_sha256")
                != execution_source_manifest_sha256
                or imported.get("destination_analysis_source_manifest_sha256")
                != destination_analysis_source_manifest_sha256
            ):
                raise SystemExit("RECOUNT_IMPORTED provenance differs from the count manifest")
    pinned_environment_hashes = {row["pinned_environment_sha256"] for row in count_rows}
    pinned_featurecounts_hashes = {row["pinned_featurecounts_executable_sha256"] for row in count_rows}
    pinned_samtools_hashes = {row["pinned_samtools_executable_sha256"] for row in count_rows}
    pinned_samtools_versions = {row["pinned_samtools_version"] for row in count_rows}
    if (
        any(len(values) != 1 for values in (
            pinned_environment_hashes, pinned_featurecounts_hashes,
            pinned_samtools_hashes, pinned_samtools_versions,
        ))
        or any(not re.fullmatch(r"[0-9a-f]{64}", value) for value in (
            next(iter(pinned_environment_hashes)), next(iter(pinned_featurecounts_hashes)),
            next(iter(pinned_samtools_hashes)),
        ))
    ):
        raise SystemExit("Count manifest does not pin one common recount tool environment")
    pinned_environment = next(iter(pinned_environment_hashes))
    pinned_featurecounts = next(iter(pinned_featurecounts_hashes))
    pinned_samtools = next(iter(pinned_samtools_hashes))
    pinned_samtools_version = next(iter(pinned_samtools_versions))
    task_environment_paths = [root / "counts" / dataset / "environment.txt" for dataset in DIRECT]
    task_environment_paths.extend(
        root / f"counts/GSE213621/chunks/chunk{index}/environment.txt" for index in range(4)
    )
    for path in task_environment_paths:
        if (
            sha256(path) != pinned_environment
            or environment_values(path, "featureCounts_executable_sha256") != {pinned_featurecounts}
            or environment_values(path, "samtools_executable_sha256") != {pinned_samtools}
            or environment_values(path, "samtools_version") != {pinned_samtools_version}
        ):
            raise SystemExit(f"Recount task environment differs from the common pin: {path}")
    for dataset, row in by_dataset.items():
        if (
            row["n_samples"] != str(EXPECTED_SAMPLES[dataset])
            or row["n_genes"] != "86369"
            or row["featurecounts_version"] != "2.1.1"
            or row["command_contract"] != "-p --countReadPairs -B -s 2"
            or row["gtf_sha256"] != EXPECTED_GTF_SHA256
        ):
            raise SystemExit(f"Count manifest scientific contract differs: {dataset}")
        commands = json.loads(row["exact_command_json"])
        expected_commands = 4 if dataset == "GSE213621" else 1
        if not isinstance(commands, list) or len(commands) != expected_commands:
            raise SystemExit(f"Count manifest exact-command cardinality differs: {dataset}")
        for command in commands:
            if not isinstance(command, str) or not command.startswith("# Program:featureCounts v2.1.1;"):
                raise SystemExit(f"Count manifest exact command lacks pinned featureCounts: {dataset}")
            for token in ('"-p"', '"--countReadPairs"', '"-B"', '"-s" "2"'):
                if token not in command:
                    raise SystemExit(f"Count manifest exact command lacks {token}: {dataset}")
        validate_dataset_artifacts(root, dataset, row, commands, manifest_samples[dataset])
        if row["bam_manifest_sha256"] != sha256(bam_manifest):
            raise SystemExit(f"Count manifest BAM digest differs: {dataset}")
        if row["source_manifest_sha256"] != execution_source_manifest_sha256:
            raise SystemExit(f"Count manifest execution-source digest differs: {dataset}")
        for field in ("featurecounts_executable_sha256_json", "samtools_executable_sha256_json"):
            values = json.loads(row[field])
            if not values or any(not re.fullmatch(r"[0-9a-f]{64}", value) for value in values):
                raise SystemExit(f"Count manifest lacks valid executable hashes: {dataset}/{field}")
        if not json.loads(row["samtools_version_json"]):
            raise SystemExit(f"Count manifest lacks samtools version: {dataset}")
        if (
            json.loads(row["featurecounts_executable_sha256_json"]) != [pinned_featurecounts]
            or json.loads(row["samtools_executable_sha256_json"]) != [pinned_samtools]
            or json.loads(row["samtools_version_json"]) != [pinned_samtools_version]
        ):
            raise SystemExit(f"Count manifest dataset tool provenance differs from the common pin: {dataset}")

    with (root / "manifests/count_overrides.tsv").open(newline="") as handle:
        overrides = list(csv.DictReader(handle, delimiter="\t"))
    if len(overrides) != 9 or {row["dataset"] for row in overrides} != ACTIVE:
        raise SystemExit("Count overrides do not contain exactly nine active cohorts")
    override_map = {row["dataset"]: Path(row["count_path"]).resolve(strict=True) for row in overrides}
    for dataset in AFFECTED:
        if override_map[dataset] != Path(by_dataset[dataset]["count_path"]).resolve(strict=True):
            raise SystemExit(f"Corrected override differs from count manifest: {dataset}")
    for dataset in ACTIVE - AFFECTED:
        expected = (root / f"frozen_sets/read_counts/{dataset}/gene_counts.txt").resolve(strict=True)
        if override_map[dataset] != expected:
            raise SystemExit(f"Unaffected override differs from frozen baseline: {dataset}")

    publication = root / "counts/GSE213621/merge_publication.json"
    merge_marker = root / "counts/GSE213621/MERGE_COMPLETE"
    if merge_marker.read_text().strip() != sha256(publication):
        raise SystemExit("GSE213621 merge publication marker differs")
    publication_payload = json.loads(publication.read_text())
    expected_publish_order = [
        "environment.txt", "featureCounts.log", "gene_counts.txt.summary",
        "validation.json", "gene_counts.txt",
    ]
    expected_hash_names = set(expected_publish_order)
    if (
        publication_payload.get("status") != "VALIDATED_STAGING"
        or publication_payload.get("run_id") != root.name
        or publication_payload.get("publish_order") != expected_publish_order
        or set(publication_payload.get("hashes", {})) != expected_hash_names
    ):
        raise SystemExit("GSE213621 merge publication contract differs")
    for name, expected in publication_payload["hashes"].items():
        if sha256(root / "counts/GSE213621" / name) != expected:
            raise SystemExit(f"GSE213621 post-publication hash drift: {name}")

    tool_versions = sorted({row["featurecounts_version"] for row in count_rows})
    if tool_versions != ["2.1.1"]:
        raise SystemExit(f"Unexpected featureCounts versions: {tool_versions}")
    all_files = expected_counts | manifests
    details = {
        "status": "PASS",
        "run_id": root.name,
        "included_bams": 820,
        "affected_cohorts": sorted(AFFECTED),
        "count_artifact_files": len(expected_counts),
        "authoritative_manifest_files": len(manifests),
        "bam_manifest_sha256": sha256(bam_manifest),
        "count_manifest_sha256": sha256(root / "manifests/count_manifest.tsv"),
        "count_overrides_sha256": sha256(root / "manifests/count_overrides.tsv"),
        "featurecounts_versions": tool_versions,
        "pinned_environment_sha256": pinned_environment,
        "pinned_featurecounts_executable_sha256": pinned_featurecounts,
        "pinned_samtools_executable_sha256": pinned_samtools,
        "pinned_samtools_version": pinned_samtools_version,
        "gtf_sha256": EXPECTED_GTF_SHA256,
        "artifact_origin": artifact_origin,
        "execution_source_run_id": execution_source_run_id,
        "execution_source_manifest_sha256": execution_source_manifest_sha256,
        "destination_analysis_source_manifest_sha256": destination_analysis_source_manifest_sha256,
    }
    return all_files, details


def parse_digest_manifest(root: Path, path: Path) -> set[Path]:
    paths: set[Path] = set()
    for line in path.read_text().splitlines():
        fields = line.split("\t", 1)
        if len(fields) != 2 or not re.fullmatch(r"[0-9a-f]{64}", fields[0]):
            raise SystemExit("Malformed recount artifact digest manifest")
        artifact = (root / fields[1]).resolve(strict=True)
        try:
            artifact.relative_to(root)
        except ValueError as exc:
            raise SystemExit(f"Recount digest path escapes run root: {artifact}") from exc
        if artifact in paths or sha256(artifact) != fields[0]:
            raise SystemExit(f"Duplicate or changed recount artifact: {artifact}")
        paths.add(artifact)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--seal", action="store_true")
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    digest_manifest = root / "manifests/recount_artifacts.sha256"
    structural = root / "comparisons/structural_validation.json"
    marker = root / "RECOUNT_COMPLETE"
    files, details = validate_semantics(root)

    if args.seal:
        legacy_temporaries = {
            digest_manifest: digest_manifest.with_suffix(".sha256.tmp"),
            structural: structural.with_suffix(".json.tmp"),
            marker: marker.with_suffix(".tmp"),
        }
        # Preflight all three publications before creating the first one so a
        # stale/symlinked legacy temp cannot leave a partially published seal.
        for path, legacy_temporary in legacy_temporaries.items():
            require_new_output_target(path, root, legacy_temporary)
        digest_text = "".join(
            f"{sha256(path)}\t{path.relative_to(root)}\n" for path in sorted(files)
        )
        atomic_write_text(
            digest_manifest,
            digest_text,
            root=root,
            legacy_temporary=legacy_temporaries[digest_manifest],
        )
        details["recount_artifacts_sha256"] = sha256(digest_manifest)
        atomic_write_text(
            structural,
            json.dumps(details, indent=2, sort_keys=True) + "\n",
            root=root,
            legacy_temporary=legacy_temporaries[structural],
        )
        marker_payload = {
            "recount_artifacts_sha256": sha256(digest_manifest),
            "structural_validation_sha256": sha256(structural),
        }
        atomic_write_text(
            marker,
            json.dumps(marker_payload, sort_keys=True) + "\n",
            root=root,
            legacy_temporary=legacy_temporaries[marker],
        )
    else:
        if not digest_manifest.is_file() or not structural.is_file() or not marker.is_file():
            raise SystemExit("Recount has not been sealed")
        marker_payload = json.loads(marker.read_text())
        if (
            marker_payload.get("recount_artifacts_sha256") != sha256(digest_manifest)
            or marker_payload.get("structural_validation_sha256") != sha256(structural)
        ):
            raise SystemExit("Recount seal digest differs")
        sealed_files = parse_digest_manifest(root, digest_manifest)
        if sealed_files != files:
            raise SystemExit("Recount sealed artifact set differs")
        report = json.loads(structural.read_text())
        if report.get("status") != "PASS" or report.get("run_id") != root.name:
            raise SystemExit("Structural validation report is not PASS for this run")
        if details["artifact_origin"] == "imported":
            imported_marker = root / "RECOUNT_IMPORTED"
            if marker_payload.get("recount_imported_sha256") != sha256(imported_marker):
                raise SystemExit("RECOUNT_COMPLETE is not cryptographically linked to RECOUNT_IMPORTED")
    print(f"PASS: verified {len(files)} immutable recount artifacts for {root.name}")


if __name__ == "__main__":
    main()
