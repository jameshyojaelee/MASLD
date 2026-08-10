#!/usr/bin/env python3
"""Fail-closed keyed merger for explicit per-sample featureCounts outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from array import array
from pathlib import Path


ANNOTATION = ("Geneid", "Chr", "Start", "End", "Strand", "Length")
DEFAULT_GENES = 86_369
COUNT_SUFFIX = ".counts.txt"


def clean_sample(value: str) -> str:
    name = Path(value).name
    suffix = ".Aligned.sortedByCoord.out.bam"
    return name[: -len(suffix)] if name.endswith(suffix) else name


def fsync_text(path: Path, lines) -> None:
    with path.open("w", newline="") as handle:
        for line in lines:
            handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def temporary_path(directory: Path, prefix: str) -> Path:
    descriptor, name = tempfile.mkstemp(dir=directory, prefix=prefix, suffix=".tmp")
    os.close(descriptor)
    return Path(name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--per-sample-dir", required=True, type=Path)
    parser.add_argument("--sample-id", action="append", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--expected-genes", type=int, default=DEFAULT_GENES)
    args = parser.parse_args()

    per_dir = args.per_sample_dir.resolve(strict=True)
    samples = args.sample_id
    if len(samples) != len(set(samples)) or any(not sample for sample in samples):
        raise SystemExit("Expected sample IDs must be unique and nonempty")
    output = args.output.resolve(strict=False)
    summary_output = args.summary.resolve(strict=False)
    if output.parent != summary_output.parent or not output.parent.is_dir():
        raise SystemExit("Merged count and summary outputs must share an existing directory")
    provenance_path = Path(f"{output}.merge_provenance.json")
    completion_path = Path(f"{output}.merge_complete")
    log_output = args.log.resolve(strict=False) if args.log else None
    if log_output and log_output.parent != output.parent:
        raise SystemExit("Transactional merge log must share the merged-output filesystem/directory")
    destinations = [output, summary_output, provenance_path, completion_path]
    if log_output:
        destinations.append(log_output)
    if any(path.exists() or path.is_symlink() for path in destinations):
        raise SystemExit("Refusing existing merged output")

    expected_counts = [per_dir / f"{sample}{COUNT_SUFFIX}" for sample in samples]
    expected_summaries = [Path(f"{path}.summary") for path in expected_counts]
    for path in expected_counts + expected_summaries:
        if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
            raise SystemExit(f"Missing, empty, or symlinked expected input: {path}")
    observed_counts = set(per_dir.glob(f"*{COUNT_SUFFIX}"))
    observed_summaries = set(per_dir.glob(f"*{COUNT_SUFFIX}.summary"))
    if observed_counts != set(expected_counts) or observed_summaries != set(expected_summaries):
        raise SystemExit("Per-sample directory contains missing or extra count/summary files")

    annotations: list[tuple[str, ...]] = []
    gene_index: dict[str, int] = {}
    count_columns: list[array] = []
    commands: list[str] = []
    paired_modes: list[bool] = []
    summary_statuses: list[str] | None = None
    summary_columns: list[list[int]] = []

    for file_index, (sample, count_path, summary_path) in enumerate(zip(samples, expected_counts, expected_summaries)):
        with count_path.open(newline="") as handle:
            command = handle.readline().rstrip("\n")
            if not command.startswith("# Program:featureCounts v2.1.1;"):
                raise SystemExit(f"Invalid featureCounts command header: {count_path}")
            paired = '"-p"' in command
            if paired and ('"--countReadPairs"' not in command or '"-B"' not in command):
                raise SystemExit(f"Paired source lacks fragment/both-end flags: {count_path}")
            commands.append(command)
            paired_modes.append(paired)
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            if tuple(header[:6]) != ANNOTATION or len(header) != 7 or clean_sample(header[6]) != sample:
                raise SystemExit(f"Count header/sample contract failure: {count_path}")
            values = array("q", [0]) * args.expected_genes
            seen: set[str] = set()
            row_count = 0
            for row in reader:
                row_count += 1
                if len(row) != 7:
                    raise SystemExit(f"Malformed count row {row_count}: {count_path}")
                gene = row[0]
                if not gene or gene in seen:
                    raise SystemExit(f"Duplicate/empty Geneid {gene}: {count_path}")
                seen.add(gene)
                annotation = tuple(row[:6])
                if file_index == 0:
                    # The first input establishes exact keyed annotations/order.
                    gene_index[gene] = row_count - 1
                    annotations.append(annotation)
                    index = row_count - 1
                else:
                    index = gene_index.get(gene, -1)
                    if index < 0 or annotations[index] != annotation:
                        raise SystemExit(f"Annotation key drift at {gene}: {count_path}")
                try:
                    value = int(row[6])
                except ValueError as exc:
                    raise SystemExit(f"Noninteger count at {gene}: {count_path}") from exc
                if value < 0:
                    raise SystemExit(f"Negative count at {gene}: {count_path}")
                values[index] = value
            if row_count != args.expected_genes or len(seen) != args.expected_genes:
                raise SystemExit(f"Expected {args.expected_genes} unique genes in {count_path}, found {row_count}")
            if len(gene_index) != args.expected_genes or set(seen) != set(gene_index):
                raise SystemExit(f"Gene set differs from first input: {count_path}")
            count_columns.append(values)

        with summary_path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            if len(header) != 2 or header[0] != "Status" or clean_sample(header[1]) != sample:
                raise SystemExit(f"Summary header/sample contract failure: {summary_path}")
            statuses: list[str] = []
            summary_values: list[int] = []
            for row in reader:
                if len(row) != 2 or not row[0] or row[0] in statuses:
                    raise SystemExit(f"Malformed or duplicate summary status: {summary_path}")
                try:
                    value = int(row[1])
                except ValueError as exc:
                    raise SystemExit(f"Noninteger summary value: {summary_path}") from exc
                if value < 0:
                    raise SystemExit(f"Negative summary value: {summary_path}")
                statuses.append(row[0]); summary_values.append(value)
            if summary_statuses is None:
                summary_statuses = statuses
            elif statuses != summary_statuses:
                raise SystemExit(f"Summary statuses/order differ: {summary_path}")
            summary_columns.append(summary_values)

    if summary_statuses is None or "Assigned" not in summary_statuses:
        raise SystemExit("Merged summaries lack the required Assigned status")
    if len(set(paired_modes)) != 1:
        raise SystemExit("Refusing to merge mixed single-end and paired-fragment count sources")

    command_digest = hashlib.sha256("\n".join(commands).encode()).hexdigest()
    count_tmp_handle = tempfile.NamedTemporaryFile(
        mode="w", newline="", delete=False, dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    summary_tmp_handle = tempfile.NamedTemporaryFile(
        mode="w", newline="", delete=False, dir=output.parent, prefix=f".{summary_output.name}.", suffix=".tmp"
    )
    count_tmp = Path(count_tmp_handle.name); summary_tmp = Path(summary_tmp_handle.name)
    try:
        with count_tmp_handle as handle:
            handle.write(
                f"# Program:featureCounts v2.1.1; exact keyed per-sample merge; "
                f"n={len(samples)}; source_commands_sha256={command_digest}\n"
            )
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow((*ANNOTATION, *samples))
            for index, annotation in enumerate(annotations):
                writer.writerow((*annotation, *(column[index] for column in count_columns)))
            handle.flush(); os.fsync(handle.fileno())
        with summary_tmp_handle as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("Status", *samples))
            for index, status in enumerate(summary_statuses or []):
                writer.writerow((status, *(column[index] for column in summary_columns)))
            handle.flush(); os.fsync(handle.fileno())
        if any(path.exists() for path in destinations):
            raise SystemExit("Merged destination appeared during validation; refusing overwrite")
        count_hash = sha256_file(count_tmp)
        summary_hash = sha256_file(summary_tmp)
        provenance = {
            "status": "PASS",
            "n_genes": len(annotations),
            "n_samples": len(samples),
            "samples": samples,
            "all_paired": all(paired_modes),
            "source_commands_sha256": command_digest,
            "count_path": str(output),
            "summary_path": str(summary_output),
            "count_sha256": count_hash,
            "summary_sha256": summary_hash,
        }
        provenance_tmp = temporary_path(output.parent, f".{provenance_path.name}.")
        log_tmp = temporary_path(output.parent, ".merge_log.") if log_output else None
        marker_tmp = temporary_path(output.parent, ".merge_complete.")
        try:
            fsync_text(provenance_tmp, [json.dumps(provenance, indent=2, sort_keys=True) + "\n"])
            if log_tmp:
                fsync_text(log_tmp, [f"PASS: exact keyed merge of {len(samples)} samples and {len(annotations)} genes\n"])
            fsync_text(marker_tmp, [json.dumps({"count_sha256": count_hash, "summary_sha256": summary_hash}, sort_keys=True) + "\n"])
            os.replace(summary_tmp, summary_output)
            os.replace(count_tmp, output)
            if sha256_file(output) != count_hash:
                raise SystemExit("Post-publication count hash differs")
            if sha256_file(summary_output) != summary_hash:
                raise SystemExit("Post-publication summary hash differs")
            os.replace(provenance_tmp, provenance_path)
            if log_tmp and log_output:
                os.replace(log_tmp, log_output)
            # This marker is the transaction commit and therefore must be last.
            os.replace(marker_tmp, completion_path)
        finally:
            provenance_tmp.unlink(missing_ok=True)
            marker_tmp.unlink(missing_ok=True)
            if log_tmp:
                log_tmp.unlink(missing_ok=True)
    finally:
        count_tmp.unlink(missing_ok=True)
        summary_tmp.unlink(missing_ok=True)

    print(json.dumps(provenance, sort_keys=True))


if __name__ == "__main__":
    main()
