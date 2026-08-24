#!/usr/bin/env python3
"""Split raw GSE296875 fragments into donor-by-lineage files without leakage."""

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
import subprocess
import sys
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file


PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
SOURCE_PRIMARY_CONTIGS = tuple(sorted(PRIMARY_CONTIGS))
CONTIG_RANK = {
    contig: index for index, contig in enumerate(SOURCE_PRIMARY_CONTIGS)
}
MEMBERSHIP_FIELDS = (
    "well_id",
    "raw_barcode",
    "cell_id",
    "donor_id",
    "source_label",
    "lineage_id",
    "analysis_role",
    "outer_fold",
)


class FragmentSplitError(ValueError):
    """Raised when a fragment source or split output violates its contract."""


class HashingReader:
    """Hash compressed bytes while a GzipFile consumes them sequentially."""

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


class HashingWriter:
    """Hash deterministic gzip bytes as they are written."""

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


def _open_deterministic_gzip(path: Path) -> tuple[io.TextIOWrapper, HashingWriter, io.BufferedWriter]:
    raw = path.open("xb")
    hashed = HashingWriter(raw)
    compressed = gzip.GzipFile(
        filename="",
        mode="wb",
        compresslevel=1,
        fileobj=hashed,
        mtime=0,
    )
    return io.TextIOWrapper(compressed, encoding="utf-8", newline=""), hashed, raw


def _close_deterministic_gzip(
    text: io.TextIOWrapper, hashed: HashingWriter, raw: io.BufferedWriter
) -> str:
    text.close()
    hashed.flush()
    raw.close()
    return hashed.digest.hexdigest()


def _read_membership(path: Path, expected_well: str) -> dict[str, dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MEMBERSHIP_FIELDS:
            raise FragmentSplitError("membership fields differ")
        rows = [dict(row) for row in reader]
    if not rows or {row["well_id"] for row in rows} != {expected_well}:
        raise FragmentSplitError(f"membership well differs for {expected_well}")
    by_barcode = {row["raw_barcode"]: row for row in rows}
    if len(by_barcode) != len(rows):
        raise FragmentSplitError(f"membership barcodes are duplicated in {expected_well}")
    return by_barcode


def _input_artifact_record(
    manifest: Mapping[str, Any], relative_path: str
) -> tuple[str, int]:
    rows = manifest.get("artifacts")
    if not isinstance(rows, list):
        raise FragmentSplitError("acquisition artifact roster is invalid")
    matches = [row for row in rows if row.get("path") == relative_path]
    if len(matches) != 1:
        raise FragmentSplitError(f"fragment artifact is not unique: {relative_path}")
    return str(matches[0]["sha256"]), int(matches[0]["size_bytes"])


def _split_well(arguments: Mapping[str, str]) -> dict[str, Any]:
    well = arguments["well"]
    fragment = Path(arguments["fragment"])
    membership = _read_membership(Path(arguments["membership"]), well)
    output = Path(arguments["output"])
    output.mkdir(mode=0o750, parents=True)
    expected_size = int(arguments["expected_size"])
    expected_sha256 = arguments["expected_sha256"]
    if fragment.is_symlink() or fragment.stat().st_size != expected_size:
        raise FragmentSplitError(f"source fragment size differs for {well}")

    writers: dict[
        tuple[str, str], tuple[io.TextIOWrapper, HashingWriter, io.BufferedWriter]
    ] = {}
    stats: dict[tuple[str, str], dict[str, int]] = {}
    seen_barcodes: set[str] = set()
    source_records = 0
    source_read_support = 0
    excluded_records = 0
    excluded_read_support = 0
    nonprimary_records = 0
    header_lines = 0
    header_values: dict[str, str] = {}
    seen_data = False
    previous: tuple[int, int, int] | None = None
    raw_handle = fragment.open("rb")
    hashing_reader = HashingReader(raw_handle)
    compressed = gzip.GzipFile(fileobj=hashing_reader, mode="rb")
    text_input = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    try:
        for line_number, line in enumerate(text_input, start=1):
            if line.startswith("#"):
                if seen_data:
                    raise FragmentSplitError(
                        f"{well} header appears after fragment data"
                    )
                header_lines += 1
                payload = line[1:].strip()
                if "=" in payload:
                    key, value = payload.split("=", 1)
                    header_values[key] = value
                continue
            seen_data = True
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise FragmentSplitError(
                    f"{well} fragment row {line_number} does not have five fields"
                )
            contig, start_text, end_text, barcode, count_text = fields
            try:
                start, end, count = int(start_text), int(end_text), int(count_text)
            except ValueError as error:
                raise FragmentSplitError(
                    f"{well} fragment row {line_number} is not integer-valued"
                ) from error
            if start < 0 or end <= start or count <= 0 or not barcode:
                raise FragmentSplitError(f"{well} fragment row {line_number} is invalid")
            source_records += 1
            source_read_support += count
            rank = CONTIG_RANK.get(contig)
            if rank is None:
                nonprimary_records += 1
                continue
            coordinate = (rank, start, end)
            if previous is not None and coordinate < previous:
                raise FragmentSplitError(f"{well} fragment source is not coordinate sorted")
            previous = coordinate
            record = membership.get(barcode)
            if record is None:
                excluded_records += 1
                excluded_read_support += count
                continue
            seen_barcodes.add(barcode)
            key = (record["donor_id"], record["lineage_id"])
            if key not in writers:
                path = output / f"donor_{key[0]}__{key[1]}.tsv.gz"
                writers[key] = _open_deterministic_gzip(path)
                stats[key] = {"records": 0, "read_support": 0}
            writers[key][0].write(
                "\t".join((contig, start_text, end_text, record["cell_id"], count_text))
                + "\n"
            )
            stats[key]["records"] += 1
            stats[key]["read_support"] += count
    finally:
        text_input.close()
        raw_handle.close()
        for text, hashed, raw in writers.values():
            if not text.closed:
                _close_deterministic_gzip(text, hashed, raw)

    if hashing_reader.digest.hexdigest() != expected_sha256:
        raise FragmentSplitError(f"source fragment SHA-256 differs for {well}")
    if (
        header_lines == 0
        or header_values.get("pipeline_name") != "cellranger-arc"
        or header_values.get("pipeline_version") != "cellranger-arc-2.0.0"
        or header_values.get("reference_version") != "2020-A"
    ):
        raise FragmentSplitError(f"source fragment header differs for {well}")
    missing_barcodes = sorted(set(membership).difference(seen_barcodes))
    if missing_barcodes:
        raise FragmentSplitError(
            f"{well} has {len(missing_barcodes)} filtered nuclei without fragments"
        )
    result = {
        "well_id": well,
        "source_records": source_records,
        "source_read_support": source_read_support,
        "excluded_records": excluded_records,
        "excluded_read_support": excluded_read_support,
        "nonprimary_records": nonprimary_records,
        "header_lines": header_lines,
        "source_pipeline": header_values["pipeline_name"],
        "source_pipeline_version": header_values["pipeline_version"],
        "source_reference_version": header_values["reference_version"],
        "selected_records": sum(row["records"] for row in stats.values()),
        "selected_read_support": sum(
            row["read_support"] for row in stats.values()
        ),
        "selected_nuclei": len(seen_barcodes),
        "groups": {
            f"{donor}\0{lineage}": value
            for (donor, lineage), value in sorted(stats.items())
        },
    }
    (output / "split_stats.json").write_text(
        json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return result


def _line_key(line: str) -> tuple[tuple[int, int, int], list[str]]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 5 or fields[0] not in CONTIG_RANK:
        raise FragmentSplitError("temporary fragment row violates primary contract")
    coordinate = (CONTIG_RANK[fields[0]], int(fields[1]), int(fields[2]))
    return coordinate, fields


def _merge_group(arguments: Mapping[str, Any]) -> dict[str, Any]:
    chunks = [Path(value) for value in arguments["chunks"]]
    output = Path(arguments["output"])
    output.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
    handles = [gzip.open(path, "rt", encoding="utf-8", newline="") for path in chunks]
    previous_by_stream: list[tuple[int, int, int] | None] = [None] * len(handles)
    heap: list[tuple[tuple[int, int, int], int, str, list[str]]] = []
    for index, handle in enumerate(handles):
        line = handle.readline()
        if line:
            coordinate, fields = _line_key(line)
            heapq.heappush(heap, (coordinate, index, line, fields))
    text, hashed, raw = _open_deterministic_gzip(output)
    records = 0
    read_support = 0
    previous_global: tuple[int, int, int] | None = None
    try:
        while heap:
            coordinate, index, line, fields = heapq.heappop(heap)
            if previous_global is not None and coordinate < previous_global:
                raise FragmentSplitError("merged fragment output is not sorted")
            previous_global = coordinate
            previous_by_stream[index] = coordinate
            text.write(line)
            records += 1
            read_support += int(fields[4])
            next_line = handles[index].readline()
            if next_line:
                next_coordinate, next_fields = _line_key(next_line)
                if next_coordinate < coordinate:
                    raise FragmentSplitError("temporary fragment chunk is not sorted")
                heapq.heappush(
                    heap, (next_coordinate, index, next_line, next_fields)
                )
        output_sha256 = _close_deterministic_gzip(text, hashed, raw)
    finally:
        for handle in handles:
            handle.close()
        if not text.closed:
            text.close()
            raw.close()
    if records != int(arguments["expected_records"]) or read_support != int(
        arguments["expected_read_support"]
    ):
        raise FragmentSplitError("merged group record or read-support totals differ")
    return {
        "donor_id": arguments["donor_id"],
        "outer_fold": int(arguments["outer_fold"]),
        "lineage_id": arguments["lineage_id"],
        "analysis_role": arguments["analysis_role"],
        "nuclei": int(arguments["nuclei"]),
        "wells": arguments["wells"],
        "records": records,
        "read_support": read_support,
        "path": arguments["relative_output"],
        "size_bytes": output.stat().st_size,
        "sha256": output_sha256,
    }


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def split_and_merge(
    *,
    acquisition: Path,
    acquisition_artifacts_sha256: str,
    membership: Path,
    membership_artifacts_sha256: str,
    output: Path,
    split_workers: int,
    merge_workers: int,
    source_paths: Sequence[Path] = (),
) -> str:
    acquisition = acquisition.resolve(strict=True)
    membership = membership.resolve(strict=True)
    if output.exists() or split_workers != 8 or not 1 <= merge_workers <= 16:
        raise FragmentSplitError("invalid output or worker contract")
    if sha256_file(acquisition / "ARTIFACTS.json") != acquisition_artifacts_sha256:
        raise FragmentSplitError("fragment acquisition ARTIFACTS SHA-256 differs")
    if sha256_file(membership / "ARTIFACTS.json") != membership_artifacts_sha256:
        raise FragmentSplitError("membership ARTIFACTS SHA-256 differs")
    membership_manifest = verify_frozen_tree(membership)
    if membership_manifest.get("metadata", {}).get("cells") != 68_398:
        raise FragmentSplitError("membership cell census differs")
    acquisition_manifest = json.loads(
        (acquisition / "ARTIFACTS.json").read_text(encoding="utf-8")
    )
    well_census = _read_tsv(membership / "well_census.tsv")
    donor_lineage = _read_tsv(membership / "donor_lineage_census.tsv")
    output.mkdir(mode=0o750)
    chunk_root = output / "temporary_well_chunks"
    fragment_root = output / "fragments"
    validation = output / "validation"
    chunk_root.mkdir(mode=0o750)
    fragment_root.mkdir(mode=0o750)
    validation.mkdir(mode=0o750)

    split_arguments = []
    for row in well_census:
        relative = f"fragments/{row['fragment_filename']}"
        expected_sha256, expected_size = _input_artifact_record(
            acquisition_manifest, relative
        )
        split_arguments.append(
            {
                "well": row["well_id"],
                "fragment": str(acquisition / relative),
                "membership": str(
                    membership / "barcodes" / f"{row['well_id']}.tsv.gz"
                ),
                "output": str(chunk_root / row["well_id"]),
                "expected_sha256": expected_sha256,
                "expected_size": str(expected_size),
            }
        )
    with ProcessPoolExecutor(max_workers=split_workers) as executor:
        split_results = list(executor.map(_split_well, split_arguments))
    split_results.sort(key=lambda row: int(row["well_id"][4:]))
    (validation / "well_split_stats.json").write_text(
        json.dumps(split_results, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )

    merge_arguments = []
    for row in donor_lineage:
        if row["observed"] != "true":
            continue
        group_key = f"{row['donor_id']}\0{row['lineage_id']}"
        chunks = sorted(
            chunk_root.glob(
                f"well*/donor_{row['donor_id']}__{row['lineage_id']}.tsv.gz"
            ),
            key=lambda path: int(path.parent.name[4:]),
        )
        expected_records = sum(
            result["groups"].get(group_key, {}).get("records", 0)
            for result in split_results
        )
        expected_read_support = sum(
            result["groups"].get(group_key, {}).get("read_support", 0)
            for result in split_results
        )
        if not chunks or expected_records <= 0 or expected_read_support <= 0:
            raise FragmentSplitError(f"observed group lacks fragments: {group_key!r}")
        relative_output = (
            f"fragments/donor_{row['donor_id']}/{row['lineage_id']}.fragments.tsv.gz"
        )
        merge_arguments.append(
            {
                "chunks": [str(path) for path in chunks],
                "output": str(output / relative_output),
                "relative_output": relative_output,
                "expected_records": expected_records,
                "expected_read_support": expected_read_support,
                **row,
            }
        )
    with ProcessPoolExecutor(max_workers=merge_workers) as executor:
        merged = list(executor.map(_merge_group, merge_arguments))
    merged.sort(key=lambda row: (int(row["donor_id"]), row["lineage_id"]))
    _write_tsv(
        output / "fragment_manifest.tsv",
        (
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
        ),
        merged,
    )
    if len(merged) != 272:
        raise FragmentSplitError(f"expected 272 observed donor-lineage files; got {len(merged)}")
    if sum(row["selected_nuclei"] for row in split_results) != 68_398:
        raise FragmentSplitError("selected nucleus total differs")
    if sum(row["records"] for row in merged) != sum(
        row["selected_records"] for row in split_results
    ):
        raise FragmentSplitError("merged record total differs from well splits")
    if sum(row["read_support"] for row in merged) != sum(
        row["selected_read_support"] for row in split_results
    ):
        raise FragmentSplitError("merged read-support total differs")

    for path in sorted(chunk_root.rglob("*"), reverse=True):
        if path.is_file() and not path.is_symlink():
            path.unlink()
        elif path.is_dir():
            path.rmdir()
    chunk_root.rmdir()
    contract = {
        "schema_version": "masld-bench-gse296875-donor-lineage-fragments-v1",
        "dataset_id": "gse296875",
        "pairing": "same_nucleus",
        "biological_unit": "donor",
        "technical_batch": "well",
        "output_files": len(merged),
        "primary_contigs": list(PRIMARY_CONTIGS),
        "fragment_contig_order": "cellranger_arc_lexicographic_primary",
        "fragment_column_5": "readSupport_preserved_for_provenance",
        "assay_signal_unit": "one_unique_fragment_record_not_readSupport",
        "source_coordinates": "cellranger_arc_2_0_0_tn5_adjusted",
        "barcode_policy": "well_namespaced_cell_id",
        "unselected_barcode_policy": "excluded_not_zero",
        "nonprimary_contig_policy": "excluded_and_counted",
        "primary_lineages": [
            "cholangiocyte",
            "fibroblast",
            "hepatocyte",
            "macrophage",
            "t_cell",
        ],
        "secondary_lineages": ["b_cell", "endothelial_cell"],
        "outer_split": "donor_outer_seed_20260821_five_folds",
        "acquisition_artifacts_sha256": acquisition_artifacts_sha256,
        "membership_artifacts_sha256": membership_artifacts_sha256,
    }
    (output / "contract.json").write_text(
        json.dumps(contract, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    if source_paths:
        source_rows = []
        for source_path in source_paths:
            resolved = source_path.resolve(strict=True)
            if not resolved.is_file() or resolved.is_symlink():
                raise FragmentSplitError("source provenance path is unsafe")
            source_rows.append(f"{sha256_file(resolved)}  {resolved}\n")
        (output / "source.sha256").write_text(
            "".join(source_rows), encoding="utf-8"
        )
    environment = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"],
        check=True,
        capture_output=True,
        text=True,
    )
    (output / "pip-freeze.txt").write_text(environment.stdout, encoding="utf-8")
    return freeze_tree(
        output,
        {
            "artifact_class": "gse296875_donor_lineage_fragments",
            "dataset_id": "gse296875",
            "donors": 39,
            "lineages": 7,
            "output_files": len(merged),
            "acquisition_artifacts_sha256": acquisition_artifacts_sha256,
            "membership_artifacts_sha256": membership_artifacts_sha256,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--acquisition-artifacts-sha256", required=True)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--membership-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--split-workers", type=int, default=8)
    parser.add_argument("--merge-workers", type=int, default=16)
    parser.add_argument("--source-path", action="append", type=Path, default=[])
    arguments = parser.parse_args()
    digest = split_and_merge(
        acquisition=arguments.acquisition,
        acquisition_artifacts_sha256=arguments.acquisition_artifacts_sha256,
        membership=arguments.membership,
        membership_artifacts_sha256=arguments.membership_artifacts_sha256,
        output=arguments.output,
        split_workers=arguments.split_workers,
        merge_workers=arguments.merge_workers,
        source_paths=arguments.source_path,
    )
    print(arguments.output.resolve(strict=True))
    print(digest)


if __name__ == "__main__":
    main()
