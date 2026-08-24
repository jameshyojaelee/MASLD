#!/usr/bin/env python3
"""Audit GSE264667 HepG2 obs/var topology without reading expression values."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np


class MetadataAuditError(RuntimeError):
    """Raised when expected old-AnnData categorical topology differs."""


def strings(dataset: h5py.Dataset) -> list[str]:
    values = dataset.asstr()[...]
    return [str(value) for value in values.tolist()]


def categorical_counts(handle: h5py.File, column: str) -> dict[str, Any]:
    codes = np.asarray(handle[f"obs/{column}"][...], dtype=np.int64)
    categories = strings(handle[f"obs/__categories/{column}"])
    if codes.size != 145473 or codes.min(initial=-1) < -1 or codes.max(initial=-1) >= len(categories):
        raise MetadataAuditError(f"categorical codes differ for {column}")
    counts = np.bincount(codes[codes >= 0], minlength=len(categories))
    return {
        "category_count": len(categories),
        "missing_count": int(np.sum(codes < 0)),
        "counts": {category: int(count) for category, count in zip(categories, counts, strict=True)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.source.is_file() or args.output.exists():
        raise MetadataAuditError("metadata topology request differs")
    with h5py.File(args.source, mode="r", locking=True) as handle:
        if handle["X"].shape != (145473, 9624):
            raise MetadataAuditError("expression shape differs")
        categorical = {
            column: categorical_counts(handle, column)
            for column in ("gene", "gene_id", "gene_transcript", "sgID_AB", "transcript")
        }
        gem_groups = np.asarray(handle["obs/gem_group"][...], dtype=np.int64)
        gem_counts = Counter(int(value) for value in gem_groups.tolist())
        barcodes = strings(handle["obs/cell_barcode"])
        var_gene_ids = strings(handle["var/gene_id"])
        var_gene_names_codes = np.asarray(handle["var/gene_name"][...], dtype=np.int64)
        var_gene_name_categories = strings(handle["var/__categories/gene_name"])
        if len(barcodes) != 145473 or len(var_gene_ids) != 9624:
            raise MetadataAuditError("obs/var axis lengths differ")
        if var_gene_names_codes.min(initial=-1) < -1 or var_gene_names_codes.max(initial=-1) >= len(var_gene_name_categories):
            raise MetadataAuditError("var gene-name codes differ")
        var_gene_names = [
            var_gene_name_categories[code] if code >= 0 else ""
            for code in var_gene_names_codes.tolist()
        ]
    var_axis_payload = "\n".join(
        f"{gene_id}\t{gene_name}"
        for gene_id, gene_name in zip(var_gene_ids, var_gene_names, strict=True)
    ).encode("utf-8")
    result = {
        "schema_version": "masld-bench-gse264667-hepg2-metadata-topology-v1",
        "dataset_id": "gse264667",
        "cell_line": "HepG2",
        "num_cells": len(barcodes),
        "num_expression_features": len(var_gene_ids),
        "obs_categorical_columns": categorical,
        "gem_group_count": len(gem_counts),
        "gem_group_cell_counts": {str(key): value for key, value in sorted(gem_counts.items())},
        "cell_barcode_unique_count": len(set(barcodes)),
        "cell_barcode_plus_gem_group_unique_count": len(set(zip(barcodes, gem_groups.tolist(), strict=True))),
        "var_gene_id_unique_count": len(set(var_gene_ids)),
        "var_gene_name_unique_count": len(set(var_gene_names)),
        "var_gene_id_name_axis_sha256": hashlib.sha256(var_axis_payload).hexdigest(),
        "expression_values_read": False,
        "biological_unit_count": 1,
        "gem_groups_as_biological_replicates": False,
        "allowed_use": "training_or_descriptive_single_pool_only",
        "gse313774_accessed": False,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "metadata_topology.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
