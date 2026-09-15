#!/usr/bin/env python3
"""Materialize one author-deposited OSCAR context without treating lanes as replicates."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np


class MaterializationError(RuntimeError):
    """Raised when the read-only OSCAR processed topology differs."""


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_metadata(path: Path, context: str) -> tuple[list[dict[str, str]], dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = ["", "cell", "barcode", "sgrna", "gene", "sample", "type"]
        if reader.fieldnames != expected:
            raise MaterializationError("OSCAR metadata schema differs")
        rows = [row for row in reader if row["type"] == context]
    if not rows or len({row[""] for row in rows}) != len(rows):
        raise MaterializationError("OSCAR context metadata identifiers differ")
    by_id = {
        row[""]: "control" if row["gene"] == "Non-Targeting" else row["gene"]
        for row in rows
    }
    return rows, by_id


def merge_chunks(chunks: list[Path], destination: Path) -> None:
    with destination.open("wb") as output:
        for chunk in chunks:
            if not chunk.is_file():
                raise MaterializationError(f"missing OSCAR matrix chunk: {chunk}")
            with chunk.open("rb") as source:
                shutil.copyfileobj(source, output, length=8 * 1024 * 1024)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context", choices=("DM", "EM"), required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--chunks", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not args.metadata.is_file():
        raise MaterializationError("OSCAR materialization request differs")
    expected_chunks = 6 if args.context == "DM" else 3
    if len(args.chunks) != expected_chunks:
        raise MaterializationError("OSCAR context chunk count differs")
    rows, metadata_groups = load_metadata(args.metadata, args.context)
    args.output.mkdir(parents=True, exist_ok=False)
    merged = args.output / f"OSCAR_{args.context}_expression_matrix.csv.gz"
    merge_chunks(args.chunks, merged)
    group_names = sorted(set(metadata_groups.values()), key=lambda value: (value != "control", value))
    group_to_index = {group: index for index, group in enumerate(group_names)}
    guide_counts = Counter(
        (
            row["sample"],
            "control" if row["gene"] == "Non-Targeting" else row["gene"],
            row["barcode"],
            row["sgrna"],
        )
        for row in rows
    )
    with (args.output / "guide_family_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "target_gene", "guide_id", "guide_sequence", "num_cells"])
        for key, count in sorted(guide_counts.items()):
            writer.writerow([*key, count])
    total_values = 0
    zero_values = 0
    noninteger_values = 0
    nonfinite_values = 0
    negative_values = 0
    minimum = math.inf
    maximum = -math.inf
    feature_count = 0
    feature_ids: set[str] = set()
    with gzip.open(merged, "rt", encoding="utf-8", newline="") as matrix_handle:
        header = next(csv.reader([matrix_handle.readline().rstrip("\r\n")]))
        if len(header) < 2 or len(set(header[1:])) != len(header) - 1:
            raise MaterializationError("OSCAR matrix header differs")
        matrix_cells = header[1:]
        if set(matrix_cells) != set(metadata_groups) or len(matrix_cells) != len(rows):
            raise MaterializationError("OSCAR matrix and metadata cell axes differ")
        group_indices = np.asarray(
            [group_to_index[metadata_groups[cell]] for cell in matrix_cells], dtype=np.int64
        )
        sums_path = args.output / f"{args.context}_deposited_scale_target_sums.tsv.gz"
        with gzip.open(sums_path, "wt", encoding="utf-8", newline="") as sums_handle:
            writer = csv.writer(sums_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["feature_id", *group_names])
            for line in matrix_handle:
                feature, separator, payload = line.rstrip("\r\n").partition(",")
                if not separator or not feature or feature in feature_ids or '"' in feature:
                    raise MaterializationError("OSCAR feature row differs")
                values = np.fromstring(payload, dtype=np.float64, sep=",")
                if values.size != len(matrix_cells):
                    raise MaterializationError("OSCAR matrix row width differs")
                feature_ids.add(feature)
                feature_count += 1
                total_values += int(values.size)
                zero_values += int(np.count_nonzero(values == 0))
                noninteger_values += int(np.count_nonzero(values != np.floor(values)))
                nonfinite_values += int(np.count_nonzero(~np.isfinite(values)))
                negative_values += int(np.count_nonzero(values < 0))
                minimum = min(minimum, float(np.min(values)))
                maximum = max(maximum, float(np.max(values)))
                grouped = np.bincount(group_indices, weights=values, minlength=len(group_names))
                writer.writerow([feature, *(format(value, ".17g") for value in grouped)])
    sample_counts = Counter(row["sample"] for row in rows)
    target_guides = {row["barcode"] for row in rows if row["gene"] != "Non-Targeting"}
    control_guides = {row["barcode"] for row in rows if row["gene"] == "Non-Targeting"}
    summary = {
        "schema_version": "masld-bench-oscar-processed-context-materialization-v1",
        "dataset_id": "cra009621_oscar",
        "context": args.context,
        "cell_count": len(rows),
        "technical_lane_count": len(sample_counts),
        "technical_lane_cell_counts": dict(sorted(sample_counts.items())),
        "screen_pool_count": 1,
        "independent_biological_unit_count": "UNRESOLVED",
        "maximum_independent_units_admitted_for_inference": 1,
        "technical_lanes_as_biological_replicates": False,
        "feature_count": feature_count,
        "target_sum_column_count": len(group_names),
        "observed_target_gene_count": len({row["gene"] for row in rows if row["gene"] != "Non-Targeting"}),
        "observed_target_guide_count": len(target_guides),
        "observed_control_guide_count": len(control_guides),
        "merged_matrix_sha256": digest(merged),
        "total_values_scanned": total_values,
        "zero_value_count": zero_values,
        "noninteger_value_count": noninteger_values,
        "nonfinite_value_count": nonfinite_values,
        "negative_value_count": negative_values,
        "minimum": minimum,
        "maximum": maximum,
        "integer_nonnegative_finite_deposited_scale": bool(
            noninteger_values == 0 and nonfinite_values == 0 and negative_values == 0
        ),
        "raw_count_semantics_authoritatively_supported": False,
        "source_outcome_derived_effect_tables_used": False,
        "normalization_performed_by_materialization": False,
        "model_fit_performed": False,
        "gse313774_accessed": False,
    }
    (args.output / "materialization_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
