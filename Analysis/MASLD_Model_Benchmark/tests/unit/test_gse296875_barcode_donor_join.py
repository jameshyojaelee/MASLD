#!/usr/bin/env python3
"""Checks for the frozen GSE296875 well-namespaced barcode-to-donor join."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts import freeze_gse296875_barcode_donor_join as join


def _write(path: Path, rows: list[tuple[str, str, str, str]]) -> Path:
    lines = ["\t".join(join.REQUIRED_COLUMNS)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _row(well: str, bare: str, donor: str) -> tuple[str, str, str, str]:
    return (f"{well}_{bare}", donor, well, f"{well}_{bare}")


class GSE296875BarcodeDonorJoinTests(unittest.TestCase):
    def test_records_that_a_bare_barcode_recurs_across_wells(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write(
                Path(directory) / "join.tsv",
                [
                    _row("well1", "AAAC-1", "331"),
                    _row("well2", "AAAC-1", "342"),
                    _row("well2", "TTTG-1", "342"),
                ],
            )
            audit = join.audit_join(join.read_join(path))
            self.assertEqual(audit["nuclei"], 3)
            self.assertEqual(audit["unique_bare_10x_barcode"], 2)
            self.assertEqual(audit["bare_10x_barcodes_shared_across_wells"], 1)
            self.assertEqual(audit["maximum_wells_sharing_one_bare_10x_barcode"], 2)
            self.assertEqual(audit["nuclei_whose_bare_barcode_is_not_globally_unique"], 2)
            self.assertFalse(audit["bare_barcode_is_a_unique_key"])
            self.assertTrue(audit["well_namespacing_is_load_bearing"])

    def test_reports_whether_a_donor_spans_wells(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            nested = _write(
                Path(directory) / "nested.tsv",
                [_row("well1", "AAAC-1", "331"), _row("well2", "TTTG-1", "342")],
            )
            self.assertTrue(join.audit_join(join.read_join(nested))["donor_is_nested_within_well"])
            spanning = _write(
                Path(directory) / "spanning.tsv",
                [_row("well1", "AAAC-1", "331"), _row("well2", "TTTG-1", "331")],
            )
            audit = join.audit_join(join.read_join(spanning))
            self.assertFalse(audit["donor_is_nested_within_well"])
            self.assertEqual(audit["multi_well_donors"], {"331": ["well1", "well2"]})

    def test_rejects_a_barcode_that_is_not_well_namespaced(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write(
                Path(directory) / "join.tsv", [("AAAC-1", "331", "well1", "AAAC-1")]
            )
            with self.assertRaisesRegex(join.BarcodeJoinError, "well-namespaced"):
                join.audit_join(join.read_join(path))

    def test_rejects_a_barcode_whose_prefix_disagrees_with_its_well(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write(
                Path(directory) / "join.tsv",
                [("well2_AAAC-1", "331", "well1", "well2_AAAC-1")],
            )
            with self.assertRaisesRegex(join.BarcodeJoinError, "well prefix"):
                join.audit_join(join.read_join(path))

    def test_rejects_a_duplicate_cell_id(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write(
                Path(directory) / "join.tsv",
                [_row("well1", "AAAC-1", "331"), _row("well1", "AAAC-1", "342")],
            )
            with self.assertRaisesRegex(join.BarcodeJoinError, "not unique"):
                join.audit_join(join.read_join(path))

    def test_a_nucleus_without_a_donor_is_never_given_one(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write(
                Path(directory) / "join.tsv",
                [_row("well1", "AAAC-1", "331"), ("well1_TTTG-1", "", "well1", "well1_TTTG-1")],
            )
            with self.assertRaisesRegex(join.BarcodeJoinError, "no donor"):
                join.audit_join(join.read_join(path))

    def test_rejects_unexpected_columns(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "join.tsv"
            path.write_text("cell_id\tdonor_id\n1\t2\n", encoding="utf-8")
            with self.assertRaisesRegex(join.BarcodeJoinError, "columns differ"):
                join.read_join(path)


if __name__ == "__main__":
    unittest.main()
