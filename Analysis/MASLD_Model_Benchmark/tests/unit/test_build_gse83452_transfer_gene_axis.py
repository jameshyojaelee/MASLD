from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts.build_gse83452_transfer_gene_axis import (
    TransferGeneAxisError,
    build_axis,
    read_external_gene_ids,
)


def write_source_axis(path: Path, gene_ids: list[str]) -> None:
    lines = ["feature_index\tstable_gene_id"]
    lines.extend(f"{index}\t{gene_id}" for index, gene_id in enumerate(gene_ids))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_external_matrix(
    path: Path, gene_ids: list[str], *, columns: int = 3, key: str = "ensembl_gene_id"
) -> None:
    header = [key] + [f"GSM{2203254 + index}" for index in range(columns)]
    lines = ["\t".join(header)]
    for offset, gene_id in enumerate(gene_ids):
        lines.append(
            "\t".join([gene_id] + [f"{7.5 + offset + index:.6f}" for index in range(columns)])
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class ExternalReaderTests(unittest.TestCase):
    def test_only_the_row_key_column_is_returned(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "matrix.tsv"
            write_external_matrix(path, ["ENSG00000000003", "ENSG00000000005"])
            ids, columns = read_external_gene_ids(path, key_field="ensembl_gene_id")
            self.assertEqual(ids, ["ENSG00000000003", "ENSG00000000005"])
            self.assertEqual(columns, 3)

    def test_a_row_with_no_payload_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "matrix.tsv"
            path.write_text(
                "ensembl_gene_id\tGSM1\nENSG00000000003\n", encoding="utf-8"
            )
            with self.assertRaises(TransferGeneAxisError):
                read_external_gene_ids(path, key_field="ensembl_gene_id")

    def test_a_foreign_row_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "matrix.tsv"
            write_external_matrix(path, ["ENSG00000000003"], key="probe_id")
            with self.assertRaises(TransferGeneAxisError):
                read_external_gene_ids(path, key_field="ensembl_gene_id")


class AxisTests(unittest.TestCase):
    def _paths(self, base: Path, source: list[str], external: list[str]):
        source_path = base / "rna_feature_axis.tsv"
        external_path = base / "matrix.tsv"
        write_source_axis(source_path, source)
        write_external_matrix(external_path, external)
        return source_path, external_path

    def test_intersection_is_sorted_and_deduplicated(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            source_path, external_path = self._paths(
                base,
                ["ENSG00000000005", "ENSG00000000003", "ENSG00000000009"],
                ["ENSG00000000003", "ENSG00000000005", "ENSG00000000007"],
            )
            shared, receipt = build_axis(
                source_axis_path=source_path,
                external_matrix_path=external_path,
                key_field="ensembl_gene_id",
                expected_source_genes=3,
                expected_external_genes=3,
                expected_external_columns=3,
            )
            self.assertEqual(shared, ["ENSG00000000003", "ENSG00000000005"])
            self.assertEqual(receipt["shared_genes"], 2)
            self.assertEqual(receipt["axis_order"], "ascending_stable_gene_id")
            self.assertFalse(receipt["external_expression_values_read"])
            self.assertFalse(receipt["external_labels_read"])
            self.assertEqual(receipt["external_series"], "GSE83452")
            self.assertEqual(receipt["source_series"], "GSE135251")

    def test_a_declared_census_that_differs_from_disk_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            source_path, external_path = self._paths(
                base, ["ENSG00000000003"], ["ENSG00000000003"]
            )
            with self.assertRaises(TransferGeneAxisError):
                build_axis(
                    source_axis_path=source_path,
                    external_matrix_path=external_path,
                    key_field="ensembl_gene_id",
                    expected_source_genes=2,
                    expected_external_genes=1,
                    expected_external_columns=3,
                )

    def test_a_versioned_identifier_is_rejected_on_both_sides(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            source_path, external_path = self._paths(
                base, ["ENSG00000000003.7"], ["ENSG00000000003"]
            )
            with self.assertRaises(TransferGeneAxisError):
                build_axis(
                    source_axis_path=source_path,
                    external_matrix_path=external_path,
                    key_field="ensembl_gene_id",
                    expected_source_genes=1,
                    expected_external_genes=1,
                    expected_external_columns=3,
                )

    def test_a_disjoint_axis_pair_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            source_path, external_path = self._paths(
                base, ["ENSG00000000003"], ["ENSG00000000009"]
            )
            with self.assertRaises(TransferGeneAxisError):
                build_axis(
                    source_axis_path=source_path,
                    external_matrix_path=external_path,
                    key_field="ensembl_gene_id",
                    expected_source_genes=1,
                    expected_external_genes=1,
                    expected_external_columns=3,
                )


if __name__ == "__main__":
    unittest.main()
