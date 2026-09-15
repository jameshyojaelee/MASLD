from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

from scripts.build_gse49541_transfer_gene_axis import (
    TransferGeneAxisError,
    build_axis,
    read_external_gene_ids,
    write_axis,
)


SOURCE_IDS = ["ENSG00000000003", "ENSG00000000005", "ENSG00000004975", "ENSG00000999999"]
EXTERNAL_IDS = ["ENSG00000004975", "ENSG00000000005", "ENSG00000000003", "ENSG00000111111"]


def write_source(path: Path, ids: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("rna_feature_index", "source_feature_index", "matrix_gene_id", "stable_gene_id"))
        for index, value in enumerate(ids):
            writer.writerow((index, index, value, value))


def write_external(path: Path, ids: list[str], columns: int = 3, key: str = "ensembl_gene_id") -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([key] + [f"GSM{index}" for index in range(columns)])
        for row, value in enumerate(ids):
            writer.writerow([value] + [f"{row + column:.17g}" for column in range(columns)])


class TransferGeneAxisTests(unittest.TestCase):
    def test_intersection_is_sorted_and_label_blind(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_source(root / "source.tsv", SOURCE_IDS)
            write_external(root / "external.tsv", EXTERNAL_IDS)
            shared, receipt = build_axis(
                source_axis_path=root / "source.tsv",
                external_matrix_path=root / "external.tsv",
                key_field="ensembl_gene_id",
                expected_source_genes=4,
                expected_external_genes=4,
                expected_external_columns=3,
            )
            self.assertEqual(
                shared,
                ["ENSG00000000003", "ENSG00000000005", "ENSG00000004975"],
            )
            self.assertEqual(shared, sorted(shared))
            self.assertIs(receipt["external_expression_values_read"], False)
            self.assertIs(receipt["external_labels_read"], False)
            self.assertEqual(receipt["shared_genes"], 3)
            self.assertEqual(receipt["external_arrays"], 3)

    def test_external_reader_returns_no_expression_payload(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_external(root / "external.tsv", EXTERNAL_IDS)
            ids, columns = read_external_gene_ids(
                root / "external.tsv", key_field="ensembl_gene_id"
            )
            self.assertEqual(columns, 3)
            self.assertEqual(set(ids), set(EXTERNAL_IDS))
            # Every returned token is a gene identifier, never a value.
            self.assertTrue(all(value.startswith("ENSG") for value in ids))

    def test_declared_census_mismatch_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_source(root / "source.tsv", SOURCE_IDS)
            write_external(root / "external.tsv", EXTERNAL_IDS)
            with self.assertRaises(TransferGeneAxisError):
                build_axis(
                    source_axis_path=root / "source.tsv",
                    external_matrix_path=root / "external.tsv",
                    key_field="ensembl_gene_id",
                    expected_source_genes=5,
                    expected_external_genes=4,
                    expected_external_columns=3,
                )

    def test_versioned_identifier_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_source(root / "source.tsv", SOURCE_IDS)
            write_external(
                root / "external.tsv",
                ["ENSG00000000003.15", "ENSG00000000005", "ENSG00000004975", "ENSG00000111111"],
            )
            with self.assertRaises(TransferGeneAxisError):
                build_axis(
                    source_axis_path=root / "source.tsv",
                    external_matrix_path=root / "external.tsv",
                    key_field="ensembl_gene_id",
                    expected_source_genes=4,
                    expected_external_genes=4,
                    expected_external_columns=3,
                )

    def test_wrong_row_key_raises(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_external(root / "external.tsv", EXTERNAL_IDS, key="probe_id")
            with self.assertRaises(TransferGeneAxisError):
                read_external_gene_ids(root / "external.tsv", key_field="ensembl_gene_id")

    def test_write_axis_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "axis.tsv"
            write_axis(target, ["ENSG00000000003"])
            with self.assertRaises(FileExistsError):
                write_axis(target, ["ENSG00000000003"])


if __name__ == "__main__":
    unittest.main()
