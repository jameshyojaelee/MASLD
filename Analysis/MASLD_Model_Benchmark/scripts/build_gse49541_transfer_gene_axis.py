#!/usr/bin/env python3
"""Intersect the GSE267145 source RNA gene axis with the GSE49541 GPL570 gene axis.

Cross-platform transfer needs one shared gene axis.  Building it reads the
*identifier* column of the external gene matrix and nothing else: no intensity
value, no participant column, and no outcome ever enters this step.  That is the
same boundary the GSE260666 activation used, and it is enforced here rather than
promised - the external reader takes only field 0 of each line and raises if a
line carries no tab-delimited payload behind it.

The emitted axis is sorted by stable gene identifier so the order is a property
of the intersection and not of either platform's deposition order.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Sequence


class TransferGeneAxisError(RuntimeError):
    """Raised when the shared gene axis would leak external values."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_source_gene_ids(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "stable_gene_id" not in reader.fieldnames:
            raise TransferGeneAxisError("source RNA feature axis lacks stable_gene_id")
        ids = [row["stable_gene_id"] for row in reader]
    if len(ids) != len(set(ids)):
        raise TransferGeneAxisError("source RNA feature axis repeats a gene")
    return ids


def read_external_gene_ids(path: Path, *, key_field: str) -> tuple[list[str], int]:
    """Read only the row key column of the external gene matrix.

    The expression payload is discarded before it is ever decoded, so no
    external value can reach a preprocessing object through this reader.
    """

    ids: list[str] = []
    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if not header or header[0] != key_field:
            raise TransferGeneAxisError(
                f"external gene matrix row key is not {key_field}"
            )
        columns = len(header) - 1
        for line in handle:
            if not line.strip():
                continue
            key, tab, _payload = line.partition("\t")
            if not tab:
                raise TransferGeneAxisError("external gene matrix row has no payload")
            ids.append(key)
    if len(ids) != len(set(ids)):
        raise TransferGeneAxisError("external gene matrix repeats a gene")
    if any("." in value for value in ids):
        raise TransferGeneAxisError("external gene identifier carries a version suffix")
    return ids, columns


def build_axis(
    *,
    source_axis_path: Path,
    external_matrix_path: Path,
    key_field: str,
    expected_source_genes: int,
    expected_external_genes: int,
    expected_external_columns: int,
) -> tuple[list[str], dict[str, object]]:
    source_ids = read_source_gene_ids(source_axis_path)
    external_ids, columns = read_external_gene_ids(
        external_matrix_path, key_field=key_field
    )
    if any("." in value for value in source_ids):
        raise TransferGeneAxisError("source gene identifier carries a version suffix")
    if (
        len(source_ids) != expected_source_genes
        or len(external_ids) != expected_external_genes
        or columns != expected_external_columns
    ):
        raise TransferGeneAxisError("declared platform census differs from disk")
    shared = sorted(set(source_ids) & set(external_ids))
    if not shared:
        raise TransferGeneAxisError("source and external gene axes do not intersect")
    receipt = {
        "schema_version": "masld-bench-gse49541-transfer-gene-axis-v1",
        "status": "pass_label_blind_shared_gene_axis",
        "source_series": "GSE267145",
        "external_series": "GSE49541",
        "external_cohort_family_id": "gse31803_gse49541_fibrosis_array",
        "source_genes": len(source_ids),
        "external_genes": len(external_ids),
        "external_arrays": columns,
        "shared_genes": len(shared),
        "axis_order": "ascending_stable_gene_id",
        "external_expression_values_read": False,
        "external_labels_read": False,
        "external_participant_metadata_read": False,
        "source_labels_read": False,
        "source_axis_sha256": sha256_file(source_axis_path),
        "external_matrix_sha256": sha256_file(external_matrix_path),
        "shared_axis_sha256": hashlib.sha256(
            "\n".join(shared).encode("utf-8")
        ).hexdigest(),
    }
    return shared, receipt


def write_axis(path: Path, gene_ids: Sequence[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("feature_index", "stable_gene_id"))
        for index, gene_id in enumerate(gene_ids):
            writer.writerow((index, gene_id))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-feature-axis", required=True, type=Path)
    parser.add_argument("--external-gene-matrix", required=True, type=Path)
    parser.add_argument("--external-key-field", default="ensembl_gene_id")
    parser.add_argument("--expected-source-genes", required=True, type=int)
    parser.add_argument("--expected-external-genes", required=True, type=int)
    parser.add_argument("--expected-external-columns", required=True, type=int)
    parser.add_argument("--output-axis", required=True, type=Path)
    parser.add_argument("--output-receipt", required=True, type=Path)
    arguments = parser.parse_args()
    shared, receipt = build_axis(
        source_axis_path=arguments.source_feature_axis,
        external_matrix_path=arguments.external_gene_matrix,
        key_field=arguments.external_key_field,
        expected_source_genes=arguments.expected_source_genes,
        expected_external_genes=arguments.expected_external_genes,
        expected_external_columns=arguments.expected_external_columns,
    )
    write_axis(arguments.output_axis, shared)
    with arguments.output_receipt.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
