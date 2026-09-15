#!/usr/bin/env python3
"""Acquire and structurally audit the fixed LS-GKM matching tracks."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any


PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
UMAP_URL = "https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bedgraph.gz"
RMSK_URL = (
    "https://www.repeatmasker.org/genomes/hg38/"
    "rmsk4.0.6_dfam2.0/hg38.fa.out.gz"
)
UMAP_SIZE = 753_183_977
RMSK_SIZE = 192_485_822
UMAP_SHA256 = "6aa2e28154086a84fef75986568207f106abb09f48e138252ef9dab9dafc5331"
RMSK_SHA256 = "c0a8fc4ce71e77cc7b4d524bec50996bba17e06a38252e2ae4c726c33ae7d6a2"
UMAP_ROWS = 126_399_648
RMSK_WIDTH_COUNTS = {14: 5, 15: 5_922_173}


class MatchingTrackError(RuntimeError):
    """Raised when a remote source or its assembly surface differs."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path, header_path: Path) -> None:
    if destination.exists() or header_path.exists():
        raise MatchingTrackError("refusing to overwrite a matching-track acquisition")
    command = [
        "curl",
        "--fail",
        "--location",
        "--proto",
        "=https",
        "--tlsv1.2",
        "--retry",
        "5",
        "--connect-timeout",
        "30",
        "--dump-header",
        str(header_path),
        "--output",
        str(destination),
        url,
    ]
    completed = subprocess.run(command, stdin=subprocess.DEVNULL, check=False, shell=False)
    if completed.returncode != 0:
        raise MatchingTrackError(f"curl failed for {url}")


def load_chrom_sizes(path: Path) -> dict[str, int]:
    expected: dict[str, int] = {}
    with path.open("r", encoding="ascii") as handle:
        for line in handle:
            contig, size = line.rstrip("\n").split("\t")
            expected[contig] = int(size)
    if tuple(expected) != PRIMARY_CONTIGS:
        raise MatchingTrackError("project primary-contig roster differs")
    return expected


def validate_rmsk(path: Path, chrom_sizes: Path) -> dict[str, Any]:
    expected = load_chrom_sizes(chrom_sizes)
    counts = {contig: 0 for contig in PRIMARY_CONTIGS}
    previous_start: dict[str, int] = {}
    out_of_order_transitions = {contig: 0 for contig in PRIMARY_CONTIGS}
    width_counts: dict[int, int] = {}
    rows = 0
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields or not fields[0].isdigit():
                continue
            if len(fields) not in RMSK_WIDTH_COUNTS:
                raise MatchingTrackError(f"RepeatMasker row {line_number} has wrong width")
            width_counts[len(fields)] = width_counts.get(len(fields), 0) + 1
            contig = fields[4]
            try:
                start1, end1 = int(fields[5]), int(fields[6])
            except ValueError as error:
                raise MatchingTrackError("RepeatMasker coordinates are not integers") from error
            if start1 < 1 or end1 < start1:
                raise MatchingTrackError("RepeatMasker interval geometry differs")
            rows += 1
            if contig not in counts:
                continue
            if end1 > expected[contig]:
                raise MatchingTrackError(f"RepeatMasker interval exceeds {contig}")
            prior = previous_start.get(contig)
            if prior is not None and start1 < prior:
                out_of_order_transitions[contig] += 1
            previous_start[contig] = start1
            counts[contig] += 1
    if (
        rows != sum(RMSK_WIDTH_COUNTS.values())
        or width_counts != RMSK_WIDTH_COUNTS
        or any(value == 0 for value in counts.values())
    ):
        raise MatchingTrackError("RepeatMasker census is incomplete")
    return {
        "rows": rows,
        "primary_contig_rows": counts,
        "row_width_counts": width_counts,
        "out_of_order_transitions_by_primary_contig": out_of_order_transitions,
        "matching_interval_order": "sorted_then_unioned_by_the_input_builder",
        "source_coordinate_system": "one_based_inclusive_RepeatMasker_out",
        "matching_coordinate_conversion": "start0=query_begin1_minus_1;end0=query_end1",
        "upstream": "RepeatMasker_open_4.0.6_Dfam_2.0_hg38",
    }


def convert_validate_umap(path: Path, output: Path, chrom_sizes: Path) -> dict[str, Any]:
    import pyBigWig

    expected = load_chrom_sizes(chrom_sizes)
    if output.exists():
        raise MatchingTrackError("refusing to overwrite derived Umap bigWig")
    counts = {contig: 0 for contig in PRIMARY_CONTIGS}
    covered_bases = {contig: 0 for contig in PRIMARY_CONTIGS}
    previous_end = {contig: 0 for contig in PRIMARY_CONTIGS}
    contig_order = {contig: index for index, contig in enumerate(PRIMARY_CONTIGS)}
    previous_contig_index = -1
    rows = 0
    bigwig = pyBigWig.open(str(output), "w")
    bigwig.addHeader(list(expected.items()), maxZooms=10)
    batch_contigs: list[str] = []
    batch_starts: list[int] = []
    batch_ends: list[int] = []
    batch_values: list[float] = []

    def flush() -> None:
        if batch_contigs:
            bigwig.addEntries(
                batch_contigs,
                batch_starts,
                ends=batch_ends,
                values=batch_values,
            )
            batch_contigs.clear()
            batch_starts.clear()
            batch_ends.clear()
            batch_values.clear()

    try:
        with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith(("#", "browser", "track")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 4 or fields[0] not in expected:
                    raise MatchingTrackError(f"Umap bedGraph row {line_number} differs")
                contig = fields[0]
                try:
                    start, end, value = int(fields[1]), int(fields[2]), float(fields[3])
                except ValueError as error:
                    raise MatchingTrackError("Umap bedGraph value is malformed") from error
                current_contig_index = contig_order[contig]
                if (
                    start < 0
                    or end <= start
                    or end > expected[contig]
                    or start < previous_end[contig]
                    or current_contig_index < previous_contig_index
                    or not math.isfinite(value)
                    or not 0.0 <= value <= 1.0
                ):
                    raise MatchingTrackError(f"Umap bedGraph geometry/value differs: {contig}")
                previous_contig_index = current_contig_index
                previous_end[contig] = end
                counts[contig] += 1
                covered_bases[contig] += end - start
                rows += 1
                batch_contigs.append(contig)
                batch_starts.append(start)
                batch_ends.append(end)
                batch_values.append(value)
                if len(batch_contigs) == 100_000:
                    flush()
        flush()
    finally:
        bigwig.close()
    observed_contigs = {contig for contig, count in counts.items() if count}
    expected_observed_contigs = set(PRIMARY_CONTIGS).difference({"chrY"})
    if rows != UMAP_ROWS or observed_contigs != expected_observed_contigs:
        raise MatchingTrackError("Umap bedGraph census is incomplete")
    reopened = pyBigWig.open(str(output))
    try:
        if reopened.chroms() != expected:
            raise MatchingTrackError("derived Umap bigWig contig sizes differ")
    finally:
        reopened.close()
    return {
        "primary_contig_bounds_match_GRCh38p14": True,
        "rows": rows,
        "primary_contig_rows": counts,
        "source_observed_contigs": sorted(observed_contigs),
        "source_absent_primary_contigs": ["chrY"],
        "absent_primary_contig_semantics": "zero_mappability_and_excluded_by_the_0.80_training_window_gate",
        "primary_contig_nonzero_covered_bases": covered_bases,
        "missing_interval_value": 0.0,
        "missing_value_reason": "absent_bedGraph_intervals_are_zero;the_2018_source_also_contains_explicit_zero_rows",
        "coordinate_system": "zero_based_half_open_bedGraph",
        "measure": "Umap_v1.1.0_100mer_multi_read_mappability",
        "derived_bigwig_sha256": file_sha256(output),
        "derived_bigwig_size_bytes": output.stat().st_size,
        "derived_bigwig_writer": f"pyBigWig_{pyBigWig.__version__}",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--umap-url", default=UMAP_URL)
    parser.add_argument("--rmsk-url", default=RMSK_URL)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise MatchingTrackError("acquisition output exists")
    arguments.output.mkdir(parents=True, mode=0o750)
    source = arguments.output / "source"
    validation = arguments.output / "validation"
    source.mkdir(mode=0o750)
    validation.mkdir(mode=0o750)

    umap = source / "k100.umap.bedgraph.gz"
    umap_bigwig = source / "k100.Umap.MultiTrackMappability.bw"
    rmsk = source / "hg38.fa.out.gz"
    download(arguments.umap_url, umap, validation / "umap.headers.txt")
    download(arguments.rmsk_url, rmsk, validation / "rmsk.headers.txt")
    if (
        umap.stat().st_size != UMAP_SIZE
        or rmsk.stat().st_size != RMSK_SIZE
        or file_sha256(umap) != UMAP_SHA256
        or file_sha256(rmsk) != RMSK_SHA256
    ):
        raise MatchingTrackError("remote matching-track size differs")

    receipt = {
        "schema_version": "masld-bench-lsgkm-matching-track-acquisition-v1",
        "status": "pass_source_acquisition_and_structural_audit",
        "sources": {
            "mappability": {
                "url": arguments.umap_url,
                "size_bytes": umap.stat().st_size,
                "sha256": file_sha256(umap),
                "assembly": "UCSC_hg38_primary_contigs_compatible_with_GRCh38p14_primary_coordinates",
                "upstream_version": "Umap_v1.1.0_track_hub_2017",
                "format": "bedGraph_gzip",
            },
            "repeats": {
                "url": arguments.rmsk_url,
                "size_bytes": rmsk.stat().st_size,
                "sha256": file_sha256(rmsk),
                "assembly": "UCSC_hg38_primary_contigs_compatible_with_GRCh38p14_primary_coordinates",
                "upstream_version": "RepeatMasker_open_4.0.6_Dfam_2.0",
                "format": "RepeatMasker_out_gzip",
            },
        },
        "umap_validation": convert_validate_umap(umap, umap_bigwig, arguments.chrom_sizes),
        "repeatmasker_validation": validate_rmsk(rmsk, arguments.chrom_sizes),
        "outcomes_read": False,
        "sealed_assets_read": False,
        "biological_training_executed": False,
        "production_predictions_generated": False,
        "redistribution": "source_fetch_instructions_only_pending_source_terms_review",
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
