#!/usr/bin/env python3
"""Export frozen scGLUE smoke matrices for Seurat WNN bridge execution."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
from typing import Any


class SeuratWNNPreparationError(ValueError):
    """Raised when WNN bridge preparation inputs differ."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SeuratWNNPreparationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def load_matrix(path: Path) -> Any:
    from scipy.sparse import load_npz

    matrix = load_npz(path).tocsr()
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise SeuratWNNPreparationError(f"matrix axes differ: {path}")
    if matrix.nnz == 0 or matrix.data.min() < 0:
        raise SeuratWNNPreparationError(f"matrix values differ: {path}")
    return matrix


def export_fold(source: Path, output: Path) -> dict[str, Any]:
    from scipy.io import mmwrite

    gene_fields, genes = read_tsv(source / "genes.tsv")
    peak_fields, peaks = read_tsv(source / "peaks.tsv")
    train_fields, training = read_tsv(source / "training_rows.tsv")
    rna_fields, query_rna = read_tsv(source / "query_rna_rows.tsv")
    atac_fields, query_atac = read_tsv(source / "query_atac_rows.tsv")
    if gene_fields != ("feature_index", "gene_id") or peak_fields != (
        "feature_index",
        "peak_id",
    ):
        raise SeuratWNNPreparationError("feature schema differs")
    if train_fields != ("training_id", "rna_state", "atac_state"):
        raise SeuratWNNPreparationError("training row schema differs")
    if rna_fields != ("rna_query_id", "rna_state") or atac_fields != (
        "atac_query_id",
        "atac_state",
    ):
        raise SeuratWNNPreparationError("query row schema differs")
    if any(row["rna_state"] != "observed" for row in training + query_rna):
        raise SeuratWNNPreparationError("RNA observed state differs")
    if any(row["atac_state"] != "observed" for row in training + query_atac):
        raise SeuratWNNPreparationError("ATAC observed state differs")
    train_rna = load_matrix(source / "training_rna.npz")
    train_atac = load_matrix(source / "training_atac.npz")
    held_rna = load_matrix(source / "query_rna.npz")
    held_atac = load_matrix(source / "query_atac.npz")
    if (
        train_rna.shape != (len(training), len(genes))
        or train_atac.shape != (len(training), len(peaks))
        or held_rna.shape != (len(query_rna), len(genes))
        or held_atac.shape != (len(query_atac), len(peaks))
    ):
        raise SeuratWNNPreparationError("matrix and identifier axes differ")
    train_ids = [row["training_id"] for row in training]
    rna_ids = [row["rna_query_id"] for row in query_rna]
    atac_ids = [row["atac_query_id"] for row in query_atac]
    if (
        len(set(train_ids)) != len(train_ids)
        or len(set(rna_ids)) != len(rna_ids)
        or len(set(atac_ids)) != len(atac_ids)
        or set(rna_ids) & set(atac_ids)
    ):
        raise SeuratWNNPreparationError("query identity contract differs")
    output.mkdir(parents=True, mode=0o750)
    for name, matrix in (
        ("training_rna", train_rna),
        ("training_atac", train_atac),
        ("query_rna", held_rna),
        ("query_atac", held_atac),
    ):
        mmwrite(output / f"{name}.mtx", matrix.T, field="integer")
    for name in (
        "genes.tsv",
        "peaks.tsv",
        "training_rows.tsv",
        "query_rna_rows.tsv",
        "query_atac_rows.tsv",
    ):
        shutil.copyfile(source / name, output / name)
    return {
        "training_nuclei": len(training),
        "query_rna_nuclei": len(query_rna),
        "query_atac_nuclei": len(query_atac),
        "genes": len(genes),
        "peaks": len(peaks),
        "query_ids_disjoint": True,
        "hidden_pair_map_read": False,
    }


def prepare(inputs: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise SeuratWNNPreparationError("output exists")
    output.mkdir(parents=True, mode=0o750)
    folds = {
        str(fold): export_fold(inputs / f"fold_{fold}", output / f"fold_{fold}")
        for fold in range(5)
    }
    receipt = {
        "schema_version": "masld-bench-seurat-wnn-bridge-prepared-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "pairing_topology": "same_nucleus",
        "outer_unit": "donor",
        "folds": folds,
        "training_pairing_available": True,
        "query_modality_ids_disjoint": True,
        "query_atac_row_order_permuted": True,
        "hidden_pair_map_read": False,
        "outcomes_read": False,
        "test_outcomes_read": False,
        "rna_to_atac_prediction_executed": False,
        "sealed_rna_conditioned_atac_eligible": False,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    prepare(arguments.inputs, arguments.output)


if __name__ == "__main__":
    main()
