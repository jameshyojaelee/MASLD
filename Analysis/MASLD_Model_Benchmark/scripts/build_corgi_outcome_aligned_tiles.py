#!/usr/bin/env python3
"""Build outcome-free Corgi tiles covering fixed development cCRE windows."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
from typing import BinaryIO, Iterable, Sequence


INPUT_BP = 524_288
FLANK_BP = 65_536
OUTPUT_BP = 393_216
BIN_BP = 64
OUTPUT_BINS = 6_144
FOLDS = tuple(range(5))
CCRE_FIELDS = (
    "contig",
    "output_start",
    "output_end",
    "window_id",
    "genomic_fold",
    "window_class",
    "ccre_class",
    "ccre_id",
    "ccre_start",
    "ccre_end",
    "input_start",
    "input_end",
    "selection_hash",
)
TILE_FIELDS = (
    "tile_id",
    "genomic_fold",
    "contig",
    "input_start",
    "input_end",
    "central_start",
    "central_end",
    "windows",
    "sequence_sha256",
)
MAP_FIELDS = CCRE_FIELDS + (
    "tile_id",
    "tile_central_start",
    "tile_central_end",
    "local_start",
    "local_end",
    "first_output_bin",
    "past_last_output_bin",
)
EXCLUDED_FIELDS = CCRE_FIELDS + ("exclusion_reason",)


class CorgiTileError(ValueError):
    """Raised when a reference or fixed-window requirement differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise CorgiTileError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


class IndexedFasta:
    """Minimal read-only FASTA accessor using a samtools-compatible .fai."""

    def __init__(self, fasta: Path, fai: Path) -> None:
        self.fasta = fasta.resolve(strict=True)
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.resolve(strict=True).read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) < 5 or fields[0] in self.index:
                raise CorgiTileError("FASTA index differs")
            self.index[fields[0]] = tuple(int(value) for value in fields[1:5])
        if not self.index:
            raise CorgiTileError("FASTA index is empty")
        self.handle: BinaryIO | None = None

    def __enter__(self) -> "IndexedFasta":
        self.handle = self.fasta.open("rb")
        return self

    def __exit__(self, *_: object) -> None:
        if self.handle is not None:
            self.handle.close()
        self.handle = None

    def fetch(self, contig: str, start: int, end: int) -> str:
        if self.handle is None or contig not in self.index:
            raise CorgiTileError("FASTA accessor is not open or contig is absent")
        length, offset, line_bases, line_width = self.index[contig]
        if not 0 <= start < end <= length or line_bases < 1 or line_width < line_bases:
            raise CorgiTileError("FASTA interval differs")
        pieces: list[bytes] = []
        position = start
        while position < end:
            line_index, within = divmod(position, line_bases)
            take = min(end - position, line_bases - within)
            self.handle.seek(offset + line_index * line_width + within)
            block = self.handle.read(take)
            if len(block) != take:
                raise CorgiTileError("FASTA interval is truncated")
            pieces.append(block)
            position += take
        sequence = b"".join(pieces).decode("ascii").upper()
        if len(sequence) != end - start or set(sequence) - set("ACGTN"):
            raise CorgiTileError("FASTA sequence alphabet or length differs")
        return sequence


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _emit_fasta(handle: gzip.GzipFile, name: str, sequence: str) -> None:
    handle.write(f">{name}\n".encode("ascii"))
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80].encode("ascii") + b"\n")


def build(
    *,
    split_contract: Path,
    reference_fasta: Path,
    reference_fai: Path,
    output: Path,
    expected_windows_per_fold: int = 16_000,
) -> dict[str, object]:
    if output.exists() or expected_windows_per_fold < 1:
        raise CorgiTileError("invalid output or expected-window contract")
    split_contract = split_contract.resolve(strict=True)
    rows = read_tsv(split_contract / "ccre_evaluation_windows.tsv", CCRE_FIELDS)
    fold_counts = {fold: 0 for fold in FOLDS}
    for row in rows:
        fold = int(row["genomic_fold"])
        if fold not in fold_counts:
            raise CorgiTileError("unexpected genomic fold")
        fold_counts[fold] += 1
    if any(value != expected_windows_per_fold for value in fold_counts.values()):
        raise CorgiTileError("fixed cCRE fold census differs")

    output.mkdir(mode=0o750)
    tile_rows: list[dict[str, object]] = []
    map_rows: list[dict[str, object]] = []
    excluded_rows: list[dict[str, object]] = []
    per_fold = {fold: {"tiles": 0, "scoreable_windows": 0, "excluded_windows": 0} for fold in FOLDS}
    fasta_path = output / "tiles.fa.gz"
    with IndexedFasta(reference_fasta, reference_fai) as reference, fasta_path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as fasta_out:
            tile_number = 0
            for fold in FOLDS:
                by_contig: dict[str, list[dict[str, str]]] = {}
                for row in rows:
                    if int(row["genomic_fold"]) == fold:
                        by_contig.setdefault(row["contig"], []).append(row)
                for contig in sorted(by_contig):
                    if contig not in reference.index:
                        raise CorgiTileError("cCRE contig is absent from FASTA")
                    contig_length = reference.index[contig][0]
                    eligible: list[dict[str, str]] = []
                    for row in sorted(
                        by_contig[contig],
                        key=lambda value: (int(value["output_start"]), int(value["output_end"]), value["window_id"]),
                    ):
                        start, end = int(row["output_start"]), int(row["output_end"])
                        if start < FLANK_BP:
                            excluded_rows.append({**row, "exclusion_reason": "left_reference_boundary"})
                            per_fold[fold]["excluded_windows"] += 1
                        elif end > contig_length - FLANK_BP:
                            excluded_rows.append({**row, "exclusion_reason": "right_reference_boundary"})
                            per_fold[fold]["excluded_windows"] += 1
                        else:
                            eligible.append(row)
                    cursor = 0
                    while cursor < len(eligible):
                        first = eligible[cursor]
                        first_start = int(first["output_start"])
                        latest_central_start = contig_length - FLANK_BP - OUTPUT_BP
                        central_start = min(max(first_start, FLANK_BP), latest_central_start)
                        central_end = central_start + OUTPUT_BP
                        input_start, input_end = central_start - FLANK_BP, central_end + FLANK_BP
                        if input_end - input_start != INPUT_BP:
                            raise CorgiTileError("Corgi tile length differs")
                        stop = cursor
                        while stop < len(eligible) and int(eligible[stop]["output_end"]) <= central_end:
                            stop += 1
                        assigned = eligible[cursor:stop]
                        if not assigned:
                            raise CorgiTileError("greedy tile assigned no cCRE windows")
                        tile_id = f"corgi_fold{fold}_tile{tile_number:05d}"
                        sequence = reference.fetch(contig, input_start, input_end)
                        sequence_sha = sha256(sequence.encode("ascii")).hexdigest()
                        _emit_fasta(fasta_out, tile_id, sequence)
                        tile_rows.append(
                            {
                                "tile_id": tile_id,
                                "genomic_fold": fold,
                                "contig": contig,
                                "input_start": input_start,
                                "input_end": input_end,
                                "central_start": central_start,
                                "central_end": central_end,
                                "windows": len(assigned),
                                "sequence_sha256": sequence_sha,
                            }
                        )
                        for row in assigned:
                            local_start = int(row["output_start"]) - central_start
                            local_end = int(row["output_end"]) - central_start
                            first_bin = local_start // BIN_BP
                            past_last_bin = (local_end + BIN_BP - 1) // BIN_BP
                            if not 0 <= local_start < local_end <= OUTPUT_BP or not 0 <= first_bin < past_last_bin <= OUTPUT_BINS:
                                raise CorgiTileError("window-to-output-bin mapping differs")
                            map_rows.append(
                                {
                                    **row,
                                    "tile_id": tile_id,
                                    "tile_central_start": central_start,
                                    "tile_central_end": central_end,
                                    "local_start": local_start,
                                    "local_end": local_end,
                                    "first_output_bin": first_bin,
                                    "past_last_output_bin": past_last_bin,
                                }
                            )
                        per_fold[fold]["tiles"] += 1
                        per_fold[fold]["scoreable_windows"] += len(assigned)
                        tile_number += 1
                        cursor = stop

    if len(map_rows) + len(excluded_rows) != len(rows):
        raise CorgiTileError("window coverage is not exhaustive")
    if len({row["window_id"] for row in map_rows + excluded_rows}) != len(rows):
        raise CorgiTileError("window assignment is not one-to-one")
    _write_tsv(output / "tiles.tsv", TILE_FIELDS, tile_rows)
    _write_tsv(output / "window_map.tsv", MAP_FIELDS, map_rows)
    _write_tsv(output / "excluded_windows.tsv", EXCLUDED_FIELDS, excluded_rows)
    summary: dict[str, object] = {
        "schema_version": "masld-bench-corgi-outcome-aligned-tiles-v1",
        "status": "pass_outcome_free_tile_contract",
        "dataset_id": "gse296875",
        "input_bp": INPUT_BP,
        "output_bp": OUTPUT_BP,
        "output_bins": OUTPUT_BINS,
        "bin_bp": BIN_BP,
        "tiles": len(tile_rows),
        "scoreable_windows": len(map_rows),
        "excluded_windows": len(excluded_rows),
        "per_fold": {str(key): value for key, value in per_fold.items()},
        "outcomes_or_labels_read": False,
        "model_predictions_read": False,
        "split_contract_artifacts_sha256": digest(split_contract / "ARTIFACTS.json"),
        "reference_fasta_sha256": digest(reference_fasta),
        "reference_fai_sha256": digest(reference_fai),
        "tiles_fasta_sha256": digest(fasta_path),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-contract", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-fai", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-windows-per-fold", type=int, default=16_000)
    args = parser.parse_args()
    summary = build(
        split_contract=args.split_contract,
        reference_fasta=args.reference_fasta,
        reference_fai=args.reference_fai,
        output=args.output,
        expected_windows_per_fold=args.expected_windows_per_fold,
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
