#!/usr/bin/env python3
"""Build deduplicated Cell Ranger ARC Tn5 insertion bigWigs."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
from hashlib import sha256
import heapq
import io
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pyBigWig


PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
SOURCE_PRIMARY_CONTIGS = tuple(sorted(PRIMARY_CONTIGS))
CONTIG_RANK = {
    contig: index for index, contig in enumerate(SOURCE_PRIMARY_CONTIGS)
}
FRAGMENT_MANIFEST_FIELDS = (
    "donor_id",
    "outer_fold",
    "lineage_id",
    "analysis_role",
    "nuclei",
    "wells",
    "records",
    "read_support",
    "path",
    "size_bytes",
    "sha256",
)


class Tn5BigWigError(ValueError):
    """Raised when inputs or output signal violate the Tn5 contract."""


class HashingReader:
    """Hash compressed bytes while gzip consumes them sequentially."""

    def __init__(self, handle: io.BufferedReader) -> None:
        self.handle = handle
        self.digest = sha256()

    def read(self, size: int = -1) -> bytes:
        value = self.handle.read(size)
        self.digest.update(value)
        return value

    def tell(self) -> int:
        return self.handle.tell()

    def seek(self, offset: int, whence: int = 0) -> int:
        return self.handle.seek(offset, whence)

    def seekable(self) -> bool:
        return self.handle.seekable()

    def readable(self) -> bool:
        return True


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_chrom_sizes(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise Tn5BigWigError(f"chrom-size row {line_number} differs")
            contig, size_text = fields
            size = int(size_text)
            if size <= 0 or contig in {row[0] for row in rows}:
                raise Tn5BigWigError("chrom-size contig or length is invalid")
            rows.append((contig, size))
    if tuple(contig for contig, _ in rows) != PRIMARY_CONTIGS:
        raise Tn5BigWigError("chrom-size order differs from the primary contract")
    return rows


def read_fragment_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FRAGMENT_MANIFEST_FIELDS:
            raise Tn5BigWigError("fragment manifest fields differ")
        rows = [dict(row) for row in reader]
    if len(rows) != 272:
        raise Tn5BigWigError("fragment manifest must contain 272 observed groups")
    return rows


def _flush_contig(
    bigwig: Any,
    contig: str,
    signal: dict[int, int],
    positions: list[int],
    *,
    before: int | None = None,
    batch_size: int = 250_000,
) -> tuple[int, int]:
    entries = 0
    signal_sum = 0
    while positions and (before is None or positions[0] < before):
        batch = []
        while (
            positions
            and (before is None or positions[0] < before)
            and len(batch) < batch_size
        ):
            batch.append(heapq.heappop(positions))
        values = [signal.pop(position) for position in batch]
        if any(value <= 0 for value in values):
            raise Tn5BigWigError("Tn5 signal contains a nonpositive entry")
        bigwig.addEntries(
            [contig] * len(batch),
            batch,
            ends=[position + 1 for position in batch],
            values=[float(value) for value in values],
        )
        entries += len(batch)
        signal_sum += sum(values)
    return entries, signal_sum


def _add_signal(signal: dict[int, int], positions: list[int], position: int) -> None:
    if position not in signal:
        signal[position] = 1
        heapq.heappush(positions, position)
    else:
        signal[position] += 1


def _build_one(arguments: Mapping[str, Any]) -> dict[str, Any]:
    source = Path(arguments["source"])
    output = Path(arguments["output"])
    chrom_sizes = [(str(contig), int(size)) for contig, size in arguments["chrom_sizes"]]
    size_by_contig = dict(chrom_sizes)
    if source.is_symlink() or source.stat().st_size != int(arguments["size_bytes"]):
        raise Tn5BigWigError("fragment source size differs")
    output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    if output.exists():
        raise Tn5BigWigError("refusing to overwrite Tn5 bigWig")

    bigwig_chrom_sizes = sorted(chrom_sizes, key=lambda row: CONTIG_RANK[row[0]])
    bigwig = pyBigWig.open(str(output), "w")
    bigwig.addHeader(bigwig_chrom_sizes)
    source_records = 0
    read_support = 0
    output_entries = 0
    output_signal = 0
    current_contig: str | None = None
    current_rank = -1
    current_signal: dict[int, int] = {}
    signal_positions: list[int] = []
    max_pending_positions = 0
    previous_key: tuple[int, int, int] | None = None
    raw = source.open("rb")
    hashed = HashingReader(raw)
    compressed = gzip.GzipFile(fileobj=hashed, mode="rb")
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    try:
        for line_number, line in enumerate(text, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise Tn5BigWigError(f"fragment row {line_number} differs")
            contig, start_text, end_text, _cell_id, count_text = fields
            try:
                start, end, count = int(start_text), int(end_text), int(count_text)
            except ValueError as error:
                raise Tn5BigWigError(
                    f"fragment row {line_number} is not integer-valued"
                ) from error
            rank = CONTIG_RANK.get(contig)
            if rank is None or rank < current_rank:
                raise Tn5BigWigError("fragment contig order differs")
            key = (rank, start, end)
            if previous_key is not None and key < previous_key:
                raise Tn5BigWigError("fragment coordinate order differs")
            previous_key = key
            if current_contig is None:
                current_contig = contig
                current_rank = rank
            elif contig != current_contig:
                entries, signal_sum = _flush_contig(
                    bigwig, current_contig, current_signal, signal_positions
                )
                output_entries += entries
                output_signal += signal_sum
                current_signal = {}
                signal_positions = []
                current_contig = contig
                current_rank = rank
            else:
                entries, signal_sum = _flush_contig(
                    bigwig,
                    current_contig,
                    current_signal,
                    signal_positions,
                    before=start,
                )
                output_entries += entries
                output_signal += signal_sum
            contig_size = size_by_contig[contig]
            # Cell Ranger ARC coordinates are already adjusted by +4/-5.
            # chromEnd remains BED-exclusive, so its base index is end - 1.
            left_position = start
            right_position = end - 1
            if (
                start < 0
                or end <= start
                or count <= 0
                or left_position < 0
                or left_position >= contig_size
                or end > contig_size
                or right_position < 0
                or right_position >= contig_size
            ):
                raise Tn5BigWigError(
                    f"fragment row {line_number} has an unrepresentable cut site"
                )
            _add_signal(current_signal, signal_positions, left_position)
            _add_signal(current_signal, signal_positions, right_position)
            max_pending_positions = max(max_pending_positions, len(current_signal))
            source_records += 1
            read_support += count
        if current_contig is not None:
            entries, signal_sum = _flush_contig(
                bigwig, current_contig, current_signal, signal_positions
            )
            output_entries += entries
            output_signal += signal_sum
        if current_signal or signal_positions:
            raise Tn5BigWigError("pending Tn5 signal remains after flush")
    finally:
        text.close()
        raw.close()
        bigwig.close()

    if hashed.digest.hexdigest() != arguments["sha256"]:
        raise Tn5BigWigError("fragment source SHA-256 differs")
    if source_records != int(arguments["records"]):
        raise Tn5BigWigError("fragment record total differs")
    if read_support != int(arguments["read_support"]):
        raise Tn5BigWigError("read-support total differs")
    if output_signal != 2 * source_records:
        raise Tn5BigWigError(
            "bigWig signal does not equal two insertions per unique fragment"
        )

    check = pyBigWig.open(str(output))
    try:
        observed_signal = 0.0
        observed_entries = 0
        for contig, _size in chrom_sizes:
            intervals = check.intervals(contig) or ()
            observed_entries += len(intervals)
            observed_signal += sum(
                (end - start) * value for start, end, value in intervals
            )
    finally:
        check.close()
    if observed_entries != output_entries or observed_signal != float(output_signal):
        raise Tn5BigWigError("reopened bigWig signal differs")

    relative_output = str(arguments["relative_output"])
    return {
        "donor_id": str(arguments["donor_id"]),
        "outer_fold": int(arguments["outer_fold"]),
        "lineage_id": str(arguments["lineage_id"]),
        "analysis_role": str(arguments["analysis_role"]),
        "nuclei": int(arguments["nuclei"]),
        "unique_fragments": source_records,
        "read_support": read_support,
        "tn5_insertions": output_signal,
        "nonzero_positions": output_entries,
        "max_pending_positions": max_pending_positions,
        "path": relative_output,
        "size_bytes": output.stat().st_size,
        "sha256": sha256_file(output),
    }


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_all(
    *,
    fragments: Path,
    fragments_artifacts_sha256: str,
    chrom_sizes: Path,
    chrom_sizes_sha256: str,
    output: Path,
    workers: int,
) -> list[dict[str, Any]]:
    fragments = fragments.resolve(strict=True)
    chrom_sizes = chrom_sizes.resolve(strict=True)
    if output.exists() or not 1 <= workers <= 16:
        raise Tn5BigWigError("invalid output or worker contract")
    if sha256_file(fragments / "ARTIFACTS.json") != fragments_artifacts_sha256:
        raise Tn5BigWigError("fragment ARTIFACTS SHA-256 differs")
    if not (fragments / "COMPLETE").is_file():
        raise Tn5BigWigError("fragment tree is not complete")
    if sha256_file(chrom_sizes) != chrom_sizes_sha256:
        raise Tn5BigWigError("chrom sizes SHA-256 differs")
    sizes = read_chrom_sizes(chrom_sizes)
    manifest = read_fragment_manifest(fragments / "fragment_manifest.tsv")
    output.mkdir(mode=0o750)
    bigwig_root = output / "bigwigs"
    validation = output / "validation"
    bigwig_root.mkdir(mode=0o750)
    validation.mkdir(mode=0o750)
    arguments = []
    for row in manifest:
        relative_output = (
            f"bigwigs/donor_{row['donor_id']}/{row['lineage_id']}.tn5.bw"
        )
        arguments.append(
            {
                **row,
                "source": str(fragments / row["path"]),
                "output": str(output / relative_output),
                "relative_output": relative_output,
                "chrom_sizes": sizes,
            }
        )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(_build_one, arguments))
    results.sort(key=lambda row: (int(row["donor_id"]), row["lineage_id"]))
    fields = (
        "donor_id",
        "outer_fold",
        "lineage_id",
        "analysis_role",
        "nuclei",
        "unique_fragments",
        "read_support",
        "tn5_insertions",
        "nonzero_positions",
        "max_pending_positions",
        "path",
        "size_bytes",
        "sha256",
    )
    write_tsv(output / "bigwig_manifest.tsv", fields, results)
    contract = {
        "schema_version": "masld-bench-deduplicated-tn5-bigwig-v1",
        "dataset_id": "gse296875",
        "biological_unit": "donor",
        "lineage_unit": "author_lineage_contract_v1",
        "coordinate_system": "GRCh38.p14_0_based_half_open",
        "source_coordinate_transform": "cellranger_arc_alignment_plus_4_minus_5",
        "source_coordinates_are_already_tn5_adjusted": True,
        "left_insertion_position": "fragment_start",
        "right_insertion_position": "fragment_end_minus_1_BED_exclusive",
        "additional_shift_applied": False,
        "fragment_column_5": "readSupport_audited_but_not_used_as_signal_weight",
        "assay_signal_unit": "one_unique_deduplicated_fragment_record",
        "upstream_fragment_helper_used": False,
        "upstream_helper_exclusion_reason": "shell_reinterprets_already_adjusted_endpoints_and_hides_readSupport_provenance",
        "output_signal_contract": "two_insertions_per_unique_fragment_record",
        "primary_contigs": list(PRIMARY_CONTIGS),
        "source_fragment_contig_order": "cellranger_arc_lexicographic_primary",
        "bigwig_header_contig_order": "cellranger_arc_lexicographic_primary",
        "fragments_artifacts_sha256": fragments_artifacts_sha256,
        "chrom_sizes_sha256": chrom_sizes_sha256,
        "outputs": len(results),
    }
    (output / "contract.json").write_text(
        json.dumps(contract, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragments", type=Path, required=True)
    parser.add_argument("--fragments-artifacts-sha256", required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--chrom-sizes-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    arguments = parser.parse_args()
    results = build_all(
        fragments=arguments.fragments,
        fragments_artifacts_sha256=arguments.fragments_artifacts_sha256,
        chrom_sizes=arguments.chrom_sizes,
        chrom_sizes_sha256=arguments.chrom_sizes_sha256,
        output=arguments.output,
        workers=arguments.workers,
    )
    print(json.dumps({"outputs": len(results)}, sort_keys=True))


if __name__ == "__main__":
    main()
