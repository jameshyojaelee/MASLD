#!/usr/bin/env python3
"""Pool donor-safe GSE296875 fragments and Tn5 tracks for model training."""

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
PRIMARY_LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
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


class TrainingPseudobulkError(ValueError):
    """Raised when donor-safe pooling violates its contract."""


class HashingWriter:
    def __init__(self, handle: io.BufferedWriter) -> None:
        self.handle = handle
        self.digest = sha256()
        self.name = str(handle.name)

    def write(self, value: bytes) -> int:
        self.digest.update(value)
        return self.handle.write(value)

    def flush(self) -> None:
        self.handle.flush()

    def tell(self) -> int:
        return self.handle.tell()

    def writable(self) -> bool:
        return True


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_deterministic_gzip(
    path: Path,
) -> tuple[io.TextIOWrapper, HashingWriter, io.BufferedWriter]:
    raw = path.open("xb")
    hashed = HashingWriter(raw)
    compressed = gzip.GzipFile(
        filename="", mode="wb", compresslevel=1, fileobj=hashed, mtime=0
    )
    return io.TextIOWrapper(compressed, encoding="utf-8", newline=""), hashed, raw


def close_deterministic_gzip(
    text: io.TextIOWrapper, hashed: HashingWriter, raw: io.BufferedWriter
) -> str:
    text.close()
    hashed.flush()
    raw.close()
    return hashed.digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FRAGMENT_MANIFEST_FIELDS:
            raise TrainingPseudobulkError("fragment manifest fields differ")
        rows = [dict(row) for row in reader]
    if len(rows) != 272:
        raise TrainingPseudobulkError("fragment manifest row count differs")
    return rows


def read_chrom_sizes(path: Path) -> list[tuple[str, int]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            contig, size_text = line.rstrip("\n").split("\t")
            rows.append((contig, int(size_text)))
    if tuple(contig for contig, _size in rows) != PRIMARY_CONTIGS:
        raise TrainingPseudobulkError("chromosome sizes differ")
    return rows


def parse_fragment_line(line: str) -> tuple[tuple[int, int, int], list[str]]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 5 or fields[0] not in CONTIG_RANK:
        raise TrainingPseudobulkError("fragment row differs")
    start, end, read_support = int(fields[1]), int(fields[2]), int(fields[4])
    if start < 0 or end <= start or read_support <= 0:
        raise TrainingPseudobulkError("fragment coordinates or readSupport differ")
    return (CONTIG_RANK[fields[0]], start, end), fields


def flush_signal(
    bigwig: Any,
    contig: str,
    signal: dict[int, int],
    positions: list[int],
    before: int | None = None,
) -> tuple[int, int]:
    entries = 0
    signal_sum = 0
    while positions and (before is None or positions[0] < before):
        batch = []
        while (
            positions
            and (before is None or positions[0] < before)
            and len(batch) < 250_000
        ):
            batch.append(heapq.heappop(positions))
        values = [signal.pop(position) for position in batch]
        bigwig.addEntries(
            [contig] * len(batch),
            batch,
            ends=[position + 1 for position in batch],
            values=[float(value) for value in values],
        )
        entries += len(batch)
        signal_sum += sum(values)
    return entries, signal_sum


def add_signal(signal: dict[int, int], positions: list[int], position: int) -> None:
    if position not in signal:
        signal[position] = 1
        heapq.heappush(positions, position)
    else:
        signal[position] += 1


def pool_one(arguments: Mapping[str, Any]) -> dict[str, Any]:
    sources = [Path(value) for value in arguments["sources"]]
    fragment_output = Path(arguments["fragment_output"])
    bigwig_output = Path(arguments["bigwig_output"])
    fragment_output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    handles = [gzip.open(path, "rt", encoding="utf-8", newline="") for path in sources]
    heap: list[tuple[tuple[int, int, int], int, str, list[str]]] = []
    for index, handle in enumerate(handles):
        line = handle.readline()
        if line:
            key, fields = parse_fragment_line(line)
            heapq.heappush(heap, (key, index, line, fields))
    text, hashed, raw = open_deterministic_gzip(fragment_output)
    bigwig = pyBigWig.open(str(bigwig_output), "w")
    chrom_sizes = [tuple(value) for value in arguments["chrom_sizes"]]
    size_by_contig = dict(chrom_sizes)
    bigwig.addHeader(sorted(chrom_sizes, key=lambda row: CONTIG_RANK[row[0]]))
    records = 0
    read_support = 0
    insertions = 0
    nonzero_positions = 0
    current_contig: str | None = None
    current_signal: dict[int, int] = {}
    signal_positions: list[int] = []
    max_pending_positions = 0
    previous_key: tuple[int, int, int] | None = None
    try:
        while heap:
            key, index, line, fields = heapq.heappop(heap)
            if previous_key is not None and key < previous_key:
                raise TrainingPseudobulkError("pooled fragments are not sorted")
            previous_key = key
            contig, start_text, end_text, _cell_id, support_text = fields
            start, end, support = int(start_text), int(end_text), int(support_text)
            if start >= size_by_contig[contig] or end > size_by_contig[contig]:
                raise TrainingPseudobulkError("fragment exceeds chromosome bounds")
            if current_contig is None:
                current_contig = contig
            elif contig != current_contig:
                entries, signal_sum = flush_signal(
                    bigwig, current_contig, current_signal, signal_positions
                )
                nonzero_positions += entries
                insertions += signal_sum
                current_signal = {}
                signal_positions = []
                current_contig = contig
            else:
                entries, signal_sum = flush_signal(
                    bigwig,
                    current_contig,
                    current_signal,
                    signal_positions,
                    before=start,
                )
                nonzero_positions += entries
                insertions += signal_sum
            text.write(line)
            add_signal(current_signal, signal_positions, start)
            add_signal(current_signal, signal_positions, end - 1)
            max_pending_positions = max(max_pending_positions, len(current_signal))
            records += 1
            read_support += support
            next_line = handles[index].readline()
            if next_line:
                next_key, next_fields = parse_fragment_line(next_line)
                heapq.heappush(heap, (next_key, index, next_line, next_fields))
        if current_contig is not None:
            entries, signal_sum = flush_signal(
                bigwig, current_contig, current_signal, signal_positions
            )
            nonzero_positions += entries
            insertions += signal_sum
        if current_signal or signal_positions:
            raise TrainingPseudobulkError("pending Tn5 signal remains after flush")
        fragment_sha256 = close_deterministic_gzip(text, hashed, raw)
    finally:
        for handle in handles:
            handle.close()
        if not text.closed:
            text.close()
            raw.close()
        bigwig.close()
    if records != int(arguments["expected_records"]):
        raise TrainingPseudobulkError("pooled unique-fragment count differs")
    if read_support != int(arguments["expected_read_support"]):
        raise TrainingPseudobulkError("pooled readSupport differs")
    if insertions != 2 * records:
        raise TrainingPseudobulkError("pooled Tn5 insertion count differs")
    return {
        "donor_test_fold": int(arguments["donor_test_fold"]),
        "donor_valid_fold": int(arguments["donor_valid_fold"]),
        "donor_train_folds": arguments["donor_train_folds"],
        "lineage_id": arguments["lineage_id"],
        "donors": int(arguments["donors"]),
        "nuclei": int(arguments["nuclei"]),
        "unique_fragments": records,
        "read_support": read_support,
        "tn5_insertions": insertions,
        "nonzero_positions": nonzero_positions,
        "max_pending_positions": max_pending_positions,
        "fragment_path": arguments["relative_fragment"],
        "fragment_size_bytes": fragment_output.stat().st_size,
        "fragment_sha256": fragment_sha256,
        "bigwig_path": arguments["relative_bigwig"],
        "bigwig_size_bytes": bigwig_output.stat().st_size,
        "bigwig_sha256": sha256_file(bigwig_output),
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
    if output.exists() or not 1 <= workers <= 5:
        raise TrainingPseudobulkError("invalid output or worker contract")
    if sha256_file(fragments / "ARTIFACTS.json") != fragments_artifacts_sha256:
        raise TrainingPseudobulkError("fragment ARTIFACTS SHA-256 differs")
    if sha256_file(chrom_sizes) != chrom_sizes_sha256:
        raise TrainingPseudobulkError("chromosome-size SHA-256 differs")
    manifest = read_manifest(fragments / "fragment_manifest.tsv")
    sizes = read_chrom_sizes(chrom_sizes)
    output.mkdir(mode=0o750)
    arguments = []
    for donor_test_fold in range(5):
        donor_valid_fold = (donor_test_fold + 1) % 5
        donor_train_folds = tuple(
            value
            for value in range(5)
            if value not in {donor_test_fold, donor_valid_fold}
        )
        for lineage in PRIMARY_LINEAGES:
            rows = [
                row
                for row in manifest
                if row["lineage_id"] == lineage
                and int(row["outer_fold"]) in donor_train_folds
            ]
            donors = {row["donor_id"] for row in rows}
            if not rows or len(donors) != len(rows):
                raise TrainingPseudobulkError("training donor-lineage roster differs")
            relative_root = f"donor_test_fold_{donor_test_fold}/{lineage}"
            arguments.append(
                {
                    "sources": [str(fragments / row["path"]) for row in rows],
                    "fragment_output": str(output / f"{relative_root}.fragments.tsv.gz"),
                    "bigwig_output": str(output / f"{relative_root}.tn5.bw"),
                    "relative_fragment": f"{relative_root}.fragments.tsv.gz",
                    "relative_bigwig": f"{relative_root}.tn5.bw",
                    "chrom_sizes": sizes,
                    "donor_test_fold": donor_test_fold,
                    "donor_valid_fold": donor_valid_fold,
                    "donor_train_folds": ",".join(map(str, donor_train_folds)),
                    "lineage_id": lineage,
                    "donors": len(donors),
                    "nuclei": sum(int(row["nuclei"]) for row in rows),
                    "expected_records": sum(int(row["records"]) for row in rows),
                    "expected_read_support": sum(
                        int(row["read_support"]) for row in rows
                    ),
                }
            )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(pool_one, arguments))
    results.sort(key=lambda row: (row["donor_test_fold"], row["lineage_id"]))
    fields = tuple(results[0])
    write_tsv(output / "training_pseudobulk_manifest.tsv", fields, results)
    contract = {
        "schema_version": "masld-bench-training-fold-pseudobulk-v1",
        "dataset_id": "gse296875",
        "biological_outer_unit": "donor",
        "donor_folds": 5,
        "donor_valid_fold_rule": "test_plus_1_mod_5",
        "training_folds_per_outer_split": 3,
        "held_donor_test_used": False,
        "held_donor_validation_used": False,
        "lineages": list(PRIMARY_LINEAGES),
        "pseudobulk_weighting": "raw_unique_fragment_sum_native_lane",
        "fragment_column_5": "readSupport_preserved_not_signal_weight",
        "left_insertion_position": "fragment_start",
        "right_insertion_position": "fragment_end_minus_1_BED_exclusive",
        "source_fragment_contig_order": "cellranger_arc_lexicographic_primary",
        "bigwig_header_contig_order": "cellranger_arc_lexicographic_primary",
        "outputs": len(results),
        "fragments_artifacts_sha256": fragments_artifacts_sha256,
        "chrom_sizes_sha256": chrom_sizes_sha256,
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
    parser.add_argument("--workers", type=int, default=5)
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
