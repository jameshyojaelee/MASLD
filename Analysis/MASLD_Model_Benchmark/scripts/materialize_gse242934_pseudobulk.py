#!/usr/bin/env python3
"""Materialize single-pool GSE242934 CRISPRa raw-count pseudobulk."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path


class MaterializationError(RuntimeError):
    """Raised when GSE242934 source topology differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        return [line.rstrip("\r\n") for line in handle]


def load_guide_library(path: Path) -> dict[str, tuple[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"oligo", "sequence", "gene"}:
        raise MaterializationError("sgRNA library schema differs")
    mapping = {}
    for row in rows:
        if row["oligo"] in mapping:
            raise MaterializationError(f"duplicate sgRNA: {row['oligo']}")
        mapping[row["oligo"]] = (row["sequence"], row["gene"])
    return mapping


def load_calls(path: Path, library: dict[str, tuple[str, str]]) -> dict[str, dict[str, object]]:
    by_cell: dict[str, list[tuple[str, str]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["sgrna"] not in library:
                raise MaterializationError(f"call has unknown sgRNA: {row['sgrna']}")
            sequence, gene = library[row["sgrna"]]
            if row["barcode"] != sequence or row["gene"] != gene:
                raise MaterializationError(f"call-to-library mapping differs: {row['sgrna']}")
            by_cell[row["cell"]].append((row["sgrna"], gene))
    calls: dict[str, dict[str, object]] = {}
    for cell, values in by_cell.items():
        guides = sorted({guide for guide, unused_gene in values})
        genes = {gene for unused_guide, gene in values}
        noncontrols = sorted(gene for gene in genes if gene != "Non-Targeting")
        has_control = "Non-Targeting" in genes
        if not noncontrols:
            classification, target = "control", "Non-Targeting"
        elif has_control:
            classification, target = "mixed_control_target", "|".join(noncontrols)
        elif len(noncontrols) == 1:
            classification, target = "single_gene", noncontrols[0]
        else:
            classification, target = "multi_gene", "|".join(noncontrols)
        calls[cell] = {
            "classification": classification,
            "target_gene": target,
            "guide_family": "|".join(guides),
            "num_guides": len(guides),
        }
    return calls


def aggregate(
    barcodes_path: Path,
    features_path: Path,
    matrix_path: Path,
    calls: dict[str, dict[str, object]],
) -> tuple[list[tuple[str, str, str]], list[str], dict[str, list[int]], dict[str, int], int]:
    barcodes = lines(barcodes_path)
    feature_rows = [line.split("\t") for line in lines(features_path)]
    if len(barcodes) != len(set(barcodes)) or any(len(row) != 3 for row in feature_rows):
        raise MaterializationError("expression axes differ")
    features = [(row[0], row[1], row[2]) for row in feature_rows]
    classes: Counter[str] = Counter()
    column_groups: list[str | None] = []
    for barcode in barcodes:
        call = calls.get(barcode)
        if call is None:
            classes["unassigned"] += 1
            column_groups.append(None)
        elif call["classification"] == "control":
            classes["control"] += 1
            column_groups.append("control")
        elif call["classification"] == "single_gene":
            classes["single_gene"] += 1
            column_groups.append(str(call["target_gene"]))
        else:
            classes[str(call["classification"])] += 1
            column_groups.append(None)
    groups = sorted({group for group in column_groups if group is not None})
    if "control" not in groups:
        raise MaterializationError("no non-targeting controls")
    sums = {group: [0] * len(features) for group in groups}
    observed_nnz = 0
    with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
        shape = None
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            fields = line.split()
            if shape is None:
                shape = tuple(int(value) for value in fields)
                if shape[:2] != (len(features), len(barcodes)):
                    raise MaterializationError("matrix axes differ")
                continue
            row_index, column_index, value = (int(field) for field in fields)
            observed_nnz += 1
            group = column_groups[column_index - 1]
            if group is not None:
                sums[group][row_index - 1] += value
    if shape is None or observed_nnz != shape[2]:
        raise MaterializationError("matrix nonzero count differs")
    return features, barcodes, sums, dict(classes), observed_nnz


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--guide-library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not args.input.is_dir() or not args.guide_library.is_file():
        raise MaterializationError("GSE242934 materialization request differs")
    raw = args.input / "raw"
    library = load_guide_library(args.guide_library)
    calls = load_calls(raw / "GSM7775286_BARCODE_10x_RNA-756_RNA-757.txt.gz", library)
    features, barcodes, sums, classes, observed_nnz = aggregate(
        raw / "GSM7775285_RNA_756_barcodes.tsv.gz",
        raw / "GSM7775285_RNA_756_features.tsv.gz",
        raw / "GSM7775285_RNA_756_matrix.mtx.gz",
        calls,
    )
    args.output.mkdir(parents=True, exist_ok=False)
    groups = sorted(sums, key=lambda value: (value != "control", value))
    with gzip.open(args.output / "raw_target_pseudobulk.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature_index", "gene_id", "gene_name", "feature_type", *groups])
        for index, (gene_id, gene_name, feature_type) in enumerate(features):
            writer.writerow([
                index + 1, gene_id, gene_name, feature_type,
                *(sums[group][index] for group in groups),
            ])
    barcode_set = set(barcodes)
    family_counts: Counter[tuple[str, str, str]] = Counter(
        (
            str(call["classification"]),
            str(call["target_gene"]),
            str(call["guide_family"]),
        )
        for barcode, call in calls.items()
        if barcode in barcode_set
    )
    with (args.output / "guide_family_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["classification", "target_gene", "guide_family", "num_cells"])
        for key, count in sorted(family_counts.items()):
            writer.writerow([*key, count])
    summary = {
        "schema_version": "masld-bench-gse242934-raw-pseudobulk-v1",
        "dataset_id": "gse242934",
        "perturbation_direction": "CRISPRa_activation",
        "source_acquisition_artifact_sha256": digest(args.input / "ARTIFACTS.json"),
        "guide_library_sha256": digest(args.guide_library),
        "num_expression_features": len(features),
        "num_expression_cells": len(barcodes),
        "num_matrix_nonzero_entries": observed_nnz,
        "num_called_cells": len(calls),
        "cell_assignment_counts": classes,
        "pseudobulk_columns": groups,
        "primary_single_gene_rule": "one unique non-control target and no non-targeting guide; exact guide family retained",
        "excluded_from_single_gene_pseudobulk": ["unassigned", "mixed_control_target", "multi_gene"],
        "biological_unit_count": 1,
        "allowed_use": "training_or_descriptive_technical_pilot_only",
        "replicate_inference_allowed": False,
        "normalization_performed": False,
        "feature_selection_performed": False,
        "gse313774_accessed": False,
    }
    (args.output / "materialization_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
