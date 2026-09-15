#!/usr/bin/env python3
"""Build nested-donor, crossed-genomic scBasset inputs from training fragments.

Only nested training donors contribute cells, binary fragment-overlap labels,
or peak filtering.  Fixed cCREs from the genomic-validation fold provide the
early-stopping labels for those training-cell output tasks.  Genomic-test
cCREs receive sequence only.  Held-donor ATAC is never read or exported.
"""

from __future__ import annotations

import argparse
from bisect import bisect_left
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
SOURCE_CONTIGS = tuple(sorted(PRIMARY_CONTIGS))
CONTIG_RANK = {contig: index for index, contig in enumerate(SOURCE_CONTIGS)}
LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
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
PSEUDOBULK_FIELDS = (
    "donor_test_fold",
    "donor_valid_fold",
    "donor_train_folds",
    "lineage_id",
    "donors",
    "nuclei",
    "unique_fragments",
    "read_support",
    "tn5_insertions",
    "nonzero_positions",
    "max_pending_positions",
    "fragment_path",
    "fragment_size_bytes",
    "fragment_sha256",
    "bigwig_path",
    "bigwig_size_bytes",
    "bigwig_sha256",
)
REGION_FIELDS = (
    "region_id",
    "block_id",
    "role",
    "sequence_start",
    "sequence_end",
)
SEQUENCE_LENGTH = 1344


class ScBassetInputError(ValueError):
    """Raised when scBasset fold inputs does not meet the leakage requirements."""


def _canonical_region_id(*, name: str, role: str, contig: str, center: int) -> str:
    """Preserve the common fixed-cCRE join key for held genomic regions."""

    parts = name.split("|")
    if role in {"valid", "test"}:
        if (
            len(parts) != 3
            or parts[0] != f"fixed_{role}"
            or not parts[1].startswith("ccre_")
            or not parts[2]
        ):
            raise ScBassetInputError("fixed cCRE identifier contract differs")
        return parts[1]
    if len(parts) < 2 or parts[0] != "train_macs2":
        raise ScBassetInputError("training peak identifier contract differs")
    return f"{contig}:{center}:{name}"


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ScBassetInputError(f"{path.name} fields differ")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ScBassetInputError(f"{path.name} is empty")
    return rows


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _read_regions(path: Path, fold_path: Path) -> list[dict[str, Any]]:
    role_by_contig = {
        contig: role
        for role, contigs in json.loads(fold_path.read_text(encoding="utf-8")).items()
        for contig in contigs
    }
    if set(role_by_contig) != set(PRIMARY_CONTIGS):
        raise ScBassetInputError("genomic fold contig roster differs")
    rows = []
    seen_ids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise ScBassetInputError(
                    f"region row {line_number} does not have ten fields"
                )
            contig, start_text, end_text, name = fields[:4]
            start, end, summit = int(start_text), int(end_text), int(fields[9])
            if contig not in role_by_contig or start < 0 or end <= start or summit < 0:
                raise ScBassetInputError("region geometry differs")
            role = role_by_contig[contig]
            prefix = name.split("|", 1)[0]
            expected_prefix = {
                "train": "train_macs2",
                "valid": "fixed_valid",
                "test": "fixed_test",
            }[role]
            if prefix != expected_prefix:
                raise ScBassetInputError("region name and genomic role differ")
            center = start + summit
            sequence_start = center - SEQUENCE_LENGTH // 2
            sequence_end = sequence_start + SEQUENCE_LENGTH
            if sequence_start < 0:
                raise ScBassetInputError("scBasset sequence window crosses a contig edge")
            region_id = _canonical_region_id(
                name=name,
                role=role,
                contig=contig,
                center=center,
            )
            if region_id in seen_ids:
                raise ScBassetInputError("scBasset region identifier is duplicated")
            seen_ids.add(region_id)
            rows.append(
                {
                    "region_id": region_id,
                    "block_id": contig,
                    "role": role,
                    "contig": contig,
                    "target_start": start,
                    "target_end": end,
                    "sequence_start": sequence_start,
                    "sequence_end": sequence_end,
                }
            )
    counts = {role: sum(row["role"] == role for row in rows) for role in ("train", "valid", "test")}
    if counts["train"] < 1000 or counts["valid"] != 16000 or counts["test"] != 16000:
        raise ScBassetInputError(f"scBasset region role census differs: {counts!r}")
    return rows


def _training_roster(
    membership: Path, training_folds: Sequence[int]
) -> tuple[list[dict[str, str]], dict[str, list[dict[str, str]]]]:
    allowed = set(training_folds)
    rows = [
        row
        for row in read_tsv(membership, MEMBERSHIP_FIELDS)
        if row["analysis_role"] == "primary"
        and row["lineage_id"] in LINEAGES
        and int(row["outer_fold"]) in allowed
    ]
    by_lineage = {
        lineage: sorted(
            [row for row in rows if row["lineage_id"] == lineage],
            key=lambda row: row["cell_id"],
        )
        for lineage in LINEAGES
    }
    if any(not values for values in by_lineage.values()):
        raise ScBassetInputError("a primary training lineage is empty")
    ordered = [row for lineage in LINEAGES for row in by_lineage[lineage]]
    cell_ids = [row["cell_id"] for row in ordered]
    if len(cell_ids) != len(set(cell_ids)):
        raise ScBassetInputError("training cell identifiers are duplicated")
    return ordered, by_lineage


def _source_rows(
    pseudobulk: Path, donor_test_fold: int, training_folds: Sequence[int]
) -> dict[str, dict[str, str]]:
    rows = read_tsv(pseudobulk / "training_pseudobulk_manifest.tsv", PSEUDOBULK_FIELDS)
    selected = [
        row for row in rows if int(row["donor_test_fold"]) == donor_test_fold
    ]
    matches = {row["lineage_id"]: row for row in selected}
    if len(selected) != len(LINEAGES) or len(matches) != len(selected):
        raise ScBassetInputError("pseudobulk fold/lineage rows are not one-to-one")
    if set(matches) != set(LINEAGES):
        raise ScBassetInputError("pseudobulk lineage roster differs")
    expected_folds = ",".join(map(str, training_folds))
    for row in matches.values():
        if row["donor_train_folds"] != expected_folds:
            raise ScBassetInputError("pseudobulk nested donor folds differ")
        source = pseudobulk / row["fragment_path"]
        if source.is_symlink() or source.stat().st_size != int(row["fragment_size_bytes"]):
            raise ScBassetInputError("pseudobulk fragment identity differs")
    return matches


def _interval_index(regions: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    index = {}
    for contig in PRIMARY_CONTIGS:
        values = [
            (int(row["target_start"]), int(row["target_end"]), position)
            for position, row in enumerate(regions)
            if row["contig"] == contig and row["role"] != "test"
        ]
        values.sort()
        if values:
            index[contig] = {
                "starts": [row[0] for row in values],
                "ends": [row[1] for row in values],
                "indices": [row[2] for row in values],
                "max_width": max(row[1] - row[0] for row in values),
            }
    return index


def _build_lineage_matrix(arguments: Mapping[str, Any]) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse

    lineage = str(arguments["lineage"])
    source = Path(arguments["source"])
    output = Path(arguments["output"])
    cells = list(arguments["cells"])
    cell_index = {cell_id: index for index, cell_id in enumerate(cells)}
    interval_index = arguments["interval_index"]
    accessible: list[set[int]] = [set() for _ in cells]
    seen_cells: set[str] = set()
    records = 0
    read_support = 0
    genomic_test_fragments_parsed_and_excluded = 0
    previous: tuple[int, int, int] | None = None
    test_contigs = set(arguments["test_contigs"])
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ScBassetInputError(f"{lineage} fragment row {line_number} differs")
            contig, start_text, end_text, cell_id, support_text = fields
            start, end, support = int(start_text), int(end_text), int(support_text)
            if contig not in CONTIG_RANK or start < 0 or end <= start or support <= 0:
                raise ScBassetInputError(f"{lineage} fragment geometry differs")
            key = (CONTIG_RANK[contig], start, end)
            if previous is not None and key < previous:
                raise ScBassetInputError(f"{lineage} fragment source is not sorted")
            previous = key
            local_cell = cell_index.get(cell_id)
            if local_cell is None:
                raise ScBassetInputError(f"{lineage} source contains a nontraining cell")
            seen_cells.add(cell_id)
            if contig in test_contigs:
                genomic_test_fragments_parsed_and_excluded += 1
                records += 1
                read_support += support
                continue
            regions_on_contig = interval_index.get(contig)
            if regions_on_contig is not None:
                starts = regions_on_contig["starts"]
                ends = regions_on_contig["ends"]
                region_indices = regions_on_contig["indices"]
                maximum_width = int(regions_on_contig["max_width"])
                low = bisect_left(starts, start - maximum_width + 1)
                high = bisect_left(starts, end)
                for candidate in range(low, high):
                    if ends[candidate] > start:
                        accessible[local_cell].add(region_indices[candidate])
            records += 1
            read_support += support
    if seen_cells != set(cells):
        raise ScBassetInputError(
            f"{lineage} fragment source missed {len(set(cells) - seen_cells)} training cells"
        )
    matrix_rows = []
    matrix_columns = []
    for column, values in enumerate(accessible):
        for row in sorted(values):
            matrix_rows.append(row)
            matrix_columns.append(column)
    matrix = sparse.csr_matrix(
        (
            np.ones(len(matrix_rows), dtype=np.int8),
            (np.asarray(matrix_rows), np.asarray(matrix_columns)),
        ),
        shape=(int(arguments["n_regions"]), len(cells)),
        dtype=np.int8,
    )
    sparse.save_npz(output, matrix, compressed=True)
    summary = {
        "lineage": lineage,
        "cells": len(cells),
        "unique_fragments": records,
        "read_support": read_support,
        "binary_accessible_pairs": int(matrix.nnz),
        "genomic_test_fragments_parsed_and_excluded": genomic_test_fragments_parsed_and_excluded,
        "target_semantics": "at_least_one_deduplicated_fragment_overlap",
        "read_support_used_as_signal_weight": False,
    }
    output.with_suffix(".json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def _write_sequences(
    *, fasta: Path, regions: Sequence[Mapping[str, Any]], output: Path
) -> None:
    import h5py
    import numpy as np
    import pyfaidx

    reference = pyfaidx.Fasta(
        str(fasta), as_raw=True, sequence_always_upper=True, rebuild=False
    )
    mapping = np.full(256, 255, dtype=np.uint8)
    for index, base in enumerate(b"ACGT"):
        mapping[base] = index
    try:
        with h5py.File(output, "x") as handle:
            values = handle.create_dataset(
                "X",
                shape=(len(regions), SEQUENCE_LENGTH),
                dtype="int8",
                chunks=(min(512, len(regions)), SEQUENCE_LENGTH),
            )
            for begin in range(0, len(regions), 512):
                block = regions[begin : begin + 512]
                encoded = np.empty((len(block), SEQUENCE_LENGTH), dtype=np.int8)
                for offset, row in enumerate(block):
                    sequence = str(
                        reference[row["contig"]][
                            int(row["sequence_start"]) : int(row["sequence_end"])
                        ]
                    ).upper()
                    if len(sequence) != SEQUENCE_LENGTH or set(sequence).difference("ACGT"):
                        raise ScBassetInputError(
                            f"ambiguous or truncated sequence: {row['region_id']}"
                        )
                    code = mapping[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
                    if np.any(code > 3):
                        raise ScBassetInputError("deterministic base encoding failed")
                    encoded[offset] = code.astype(np.int8)
                values[begin : begin + len(block)] = encoded
            handle.attrs["sequence_length"] = SEQUENCE_LENGTH
            handle.attrs["base_order"] = "A,C,G,T"
            handle.attrs["ambiguous_base_policy"] = "reject"
    finally:
        reference.close()


def _remove_intermediate_lineage_matrices(output: Path) -> None:
    """Remove full-axis worker matrices before the read-only output file is frozen.

    Their genomic-test rows are structurally missing and therefore must not
    survive as all-zero rows that could be mistaken for observed negatives.
    """

    directory = output / "lineage_matrices"
    expected = {
        f"{lineage}{suffix}"
        for lineage in LINEAGES
        for suffix in (".npz", ".json")
    }
    observed = {path.name for path in directory.iterdir()}
    if observed != expected or any(path.is_symlink() for path in directory.iterdir()):
        raise ScBassetInputError("intermediate lineage-matrix inventory differs")
    for name in sorted(expected):
        (directory / name).unlink()
    directory.rmdir()


def build(
    *,
    pseudobulk: Path,
    membership: Path,
    donor_folds: Path,
    regions_bed: Path,
    fold_json: Path,
    fasta: Path,
    donor_test_fold: int,
    donor_valid_fold: int,
    output: Path,
    minimum_cell_fraction: float,
    workers: int,
) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse

    if (
        output.exists()
        or donor_test_fold not in range(5)
        or donor_valid_fold not in range(5)
        or donor_test_fold == donor_valid_fold
        or not 0 < minimum_cell_fraction <= 1
        or not 1 <= workers <= 5
    ):
        raise ScBassetInputError("output, split, filter, or worker contract is invalid")
    training_folds = tuple(
        fold for fold in range(5) if fold not in {donor_test_fold, donor_valid_fold}
    )
    regions = _read_regions(regions_bed, fold_json)
    training_cells, cells_by_lineage = _training_roster(membership, training_folds)
    source_by_lineage = _source_rows(pseudobulk, donor_test_fold, training_folds)
    output.mkdir(mode=0o750)
    (output / "lineage_matrices").mkdir(mode=0o750)
    interval_index = _interval_index(regions)
    test_contigs = sorted(
        {str(row["contig"]) for row in regions if row["role"] == "test"}
    )
    worker_arguments = []
    for lineage in LINEAGES:
        row = source_by_lineage[lineage]
        worker_arguments.append(
            {
                "lineage": lineage,
                "source": str(pseudobulk / row["fragment_path"]),
                "output": str(output / "lineage_matrices" / f"{lineage}.npz"),
                "cells": [value["cell_id"] for value in cells_by_lineage[lineage]],
                "interval_index": interval_index,
                "test_contigs": test_contigs,
                "n_regions": len(regions),
            }
        )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        worker_summaries = list(executor.map(_build_lineage_matrix, worker_arguments))
    if any(
        int(row["genomic_test_fragments_parsed_and_excluded"]) <= 0
        for row in worker_summaries
    ):
        raise ScBassetInputError("a lineage did not exercise genomic-test exclusion")
    matrices = [
        sparse.load_npz(output / "lineage_matrices" / f"{lineage}.npz")
        for lineage in LINEAGES
    ]
    full = sparse.hstack(matrices, format="csr", dtype=np.int8)
    if full.shape != (len(regions), len(training_cells)) or np.any(full.data != 1):
        raise ScBassetInputError("assembled binary matrix differs")
    train_indices = np.asarray(
        [index for index, row in enumerate(regions) if row["role"] == "train"],
        dtype=np.int64,
    )
    valid_indices = np.asarray(
        [index for index, row in enumerate(regions) if row["role"] == "valid"],
        dtype=np.int64,
    )
    test_indices = np.asarray(
        [index for index, row in enumerate(regions) if row["role"] == "test"],
        dtype=np.int64,
    )
    minimum_cells = int(math.ceil(minimum_cell_fraction * len(training_cells)))
    accessible_cells = np.asarray(full[train_indices].getnnz(axis=1)).ravel()
    keep = accessible_cells >= minimum_cells
    retained_train_indices = train_indices[keep]
    if len(retained_train_indices) < 1000:
        raise ScBassetInputError("fewer than 1,000 train-only peaks pass cell filtering")
    train_matrix = full[retained_train_indices].tocsr()
    valid_matrix = full[valid_indices].tocsr()
    sparse.save_npz(output / "m_train.npz", train_matrix, compressed=True)
    sparse.save_npz(output / "m_valid.npz", valid_matrix, compressed=True)
    _remove_intermediate_lineage_matrices(output)

    selected = {
        "train": [regions[index] for index in retained_train_indices],
        "valid": [regions[index] for index in valid_indices],
        "test": [regions[index] for index in test_indices],
    }
    for role, rows in selected.items():
        _write_sequences(fasta=fasta, regions=rows, output=output / f"{role}_seqs.h5")
        write_tsv(
            output / f"{role}_regions.tsv",
            REGION_FIELDS,
            ({field: row[field] for field in REGION_FIELDS} for row in rows),
        )
    write_tsv(
        output / "training_cells.tsv",
        ("cell_id", "donor_id", "lineage", "outer_fold"),
        (
            {
                "cell_id": row["cell_id"],
                "donor_id": row["donor_id"],
                "lineage": row["lineage_id"],
                "outer_fold": row["outer_fold"],
            }
            for row in training_cells
        ),
    )
    donor_rows = read_tsv(donor_folds, ("donor_id", "outer_fold"))
    for role, fold in (("valid", donor_valid_fold), ("test", donor_test_fold)):
        rows = [
            {
                "donor_id": row["donor_id"],
                "outer_fold": row["outer_fold"],
                "evaluation_role": role,
            }
            for row in donor_rows
            if int(row["outer_fold"]) == fold
        ]
        write_tsv(
            output / f"held_{role}_donors.tsv",
            ("donor_id", "outer_fold", "evaluation_role"),
            rows,
        )
    summary = {
        "schema_version": "masld-bench-scbasset-fold-inputs-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "split_id": f"donor{donor_test_fold}_genomic{donor_test_fold}",
        "donor_train_folds": list(training_folds),
        "donor_valid_fold": donor_valid_fold,
        "donor_test_fold": donor_test_fold,
        "training_cells": len(training_cells),
        "training_donors": len({row["donor_id"] for row in training_cells}),
        "training_lineage_cells": {
            lineage: len(cells_by_lineage[lineage]) for lineage in LINEAGES
        },
        "candidate_training_peaks": len(train_indices),
        "retained_training_peaks": len(retained_train_indices),
        "minimum_training_cell_fraction": minimum_cell_fraction,
        "minimum_training_cells": minimum_cells,
        "validation_ccres": len(valid_indices),
        "test_ccres": len(test_indices),
        "train_binary_accessible_pairs": int(train_matrix.nnz),
        "valid_binary_accessible_pairs": int(valid_matrix.nnz),
        "training_target_semantics": "binary_deduplicated_fragment_overlap",
        "read_support_used_as_signal_weight": False,
        "sequence_length": SEQUENCE_LENGTH,
        "sequence_base_order": "A,C,G,T",
        "ambiguous_sequence_policy": "reject",
        "held_donor_atac_read": False,
        "held_donor_atac_exported": False,
        "source_fragment_stream_contains_genomic_test_contigs": True,
        "genomic_test_atac_parsed_only_to_exclude": True,
        "genomic_test_atac_used_for_matrix": False,
        "genomic_test_atac_exported": False,
        "model_training_input_contains_genomic_test_atac": False,
        "intermediate_lineage_matrices_removed": True,
        "missing_evidence_encoded_as_zero": False,
        "public_tutorial_weights_used": False,
        "worker_summaries": worker_summaries,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk", type=Path, required=True)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--donor-folds", type=Path, required=True)
    parser.add_argument("--regions-bed", type=Path, required=True)
    parser.add_argument("--fold-json", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--donor-test-fold", type=int, required=True)
    parser.add_argument("--donor-valid-fold", type=int, required=True)
    parser.add_argument("--minimum-cell-fraction", type=float, default=0.05)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = build(
        pseudobulk=arguments.pseudobulk,
        membership=arguments.membership,
        donor_folds=arguments.donor_folds,
        regions_bed=arguments.regions_bed,
        fold_json=arguments.fold_json,
        fasta=arguments.fasta,
        donor_test_fold=arguments.donor_test_fold,
        donor_valid_fold=arguments.donor_valid_fold,
        output=arguments.output,
        minimum_cell_fraction=arguments.minimum_cell_fraction,
        workers=arguments.workers,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
