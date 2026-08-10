#!/usr/bin/env python3
"""Fail-closed validation for one featureCounts candidate output."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path


ANNOTATION = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
EXPECTED_GENES = 86_369


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clean_sample(value: str) -> str:
    name = Path(value).name
    suffix = ".Aligned.sortedByCoord.out.bam"
    return name[: -len(suffix)] if name.endswith(suffix) else name


def expected_manifest_samples(path: Path, dataset: str, chunk_index: int | None) -> list[str]:
    with path.open(newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["dataset"] == dataset and row["expected_status"] == "included"
        ]
    expected = [row["sample_id"] for row in rows]
    if chunk_index is not None:
        start = chunk_index * 92
        expected = expected[start : min(start + 92, len(expected))]
    return expected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--expected-samples", required=True, type=int)
    parser.add_argument("--chunk-index", type=int, choices=range(4))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    for path in (args.manifest, args.counts, args.summary, args.log, args.environment):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required file: {path}")

    log_text = args.log.read_text(errors="replace")
    required_log = (
        "Paired-end reads are included",
        "Count read pairs : yes",
    )
    missing_log = [token for token in required_log if token not in log_text]
    if missing_log:
        raise SystemExit(f"featureCounts log lacks fragment-mode evidence: {missing_log}")
    if "reads are assigned on the single-end mode" in log_text:
        raise SystemExit("featureCounts log reports single-end assignment mode")
    environment_text = args.environment.read_text(errors="replace")
    for token in (
        "featureCounts_version\tfeatureCounts v2.1.1",
        "featureCounts_executable_sha256\t",
        "samtools_version\t",
        "samtools_executable_sha256\t",
        "gtf_sha256\t73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
        "micromamba_explicit_begin",
        "micromamba_explicit_end",
    ):
        if token not in environment_text:
            raise SystemExit(f"Environment record lacks required provenance: {token}")

    with args.counts.open(newline="") as handle:
        command = handle.readline().rstrip("\n")
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:6] != ANNOTATION:
            raise SystemExit(f"Unexpected annotation header: {header[:6]}")
        samples = [clean_sample(x) for x in header[6:]]
        if len(samples) != args.expected_samples or len(set(samples)) != len(samples):
            raise SystemExit(
                f"Expected {args.expected_samples} unique samples, got {len(samples)} / {len(set(samples))}"
            )
        genes: set[str] = set()
        row_count = 0
        for row in reader:
            row_count += 1
            if len(row) != len(header):
                raise SystemExit(f"Truncated count row {row_count}: {len(row)} != {len(header)}")
            gene = row[0]
            if gene in genes:
                raise SystemExit(f"Duplicate Geneid: {gene}")
            genes.add(gene)
            try:
                values = [int(x) for x in row[6:]]
            except ValueError as exc:
                raise SystemExit(f"Non-integer count at gene {gene}") from exc
            if any(value < 0 for value in values):
                raise SystemExit(f"Negative count at gene {gene}")
    if row_count != EXPECTED_GENES:
        raise SystemExit(f"Expected {EXPECTED_GENES} genes, found {row_count}")
    if not command.startswith("# Program:featureCounts v2.1.1;"):
        raise SystemExit("featureCounts command header is not exactly v2.1.1")
    for token in (
        '"-p"',
        '"--countReadPairs"',
        '"-B"',
        '"-s" "2"',
        "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz",
    ):
        if token not in command:
            raise SystemExit(f"featureCounts command header lacks required token sequence: {token}")

    with args.summary.open(newline="") as handle:
        summary_reader = csv.reader(handle, delimiter="\t")
        summary_header = next(summary_reader)
        statuses: set[str] = set()
        for row_number, row in enumerate(summary_reader, start=2):
            if len(row) != len(summary_header) or not row[0] or row[0] in statuses:
                raise SystemExit(f"Malformed or duplicate summary status at row {row_number}")
            statuses.add(row[0])
            try:
                summary_values = [int(value) for value in row[1:]]
            except ValueError as exc:
                raise SystemExit(f"Noninteger summary value at row {row_number}") from exc
            if any(value < 0 for value in summary_values):
                raise SystemExit(f"Negative summary value at row {row_number}")
    if not statuses or "Assigned" not in statuses:
        raise SystemExit("Summary lacks required unique Assigned status")
    summary_samples = [clean_sample(x) for x in summary_header[1:]]
    if summary_samples != samples:
        raise SystemExit("Count and summary sample columns differ")
    expected = expected_manifest_samples(args.manifest, args.dataset, args.chunk_index)
    if samples != expected:
        raise SystemExit("Count sample columns do not match the frozen manifest order")

    result = {
        "dataset": args.dataset,
        "command": command,
        "n_genes": row_count,
        "n_samples": len(samples),
        "samples": samples,
        "counts_sha256": digest(args.counts),
        "summary_sha256": digest(args.summary),
        "log_sha256": digest(args.log),
        "environment_sha256": digest(args.environment),
        "featurecounts_version": "2.1.1",
        "gtf_sha256": "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
        "status": "PASS",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(args.output.suffix + ".tmp")
    tmp.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, args.output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
