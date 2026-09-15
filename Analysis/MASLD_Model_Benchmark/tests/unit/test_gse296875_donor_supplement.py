#!/usr/bin/env python3
"""Checks for the fail-closed GSE296875 donor-supplement acquisition."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import tomllib
import unittest
import zipfile

from masld_bench.artifacts import verify_frozen_tree
from scripts import acquire_gse296875_donor_supplement as supplement


ROOT = Path(__file__).parents[2]
TASKSPEC = ROOT / "config/evaluation/gse296875_histopathology.toml"

_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_DOC_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PKG_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def _workbook(path: Path, sheet_name: str, grid: list[list[str | None]]) -> None:
    """Write a minimal inline-string XLSX so tests need no third-party reader."""

    def cell(reference: str, value: str | None) -> str:
        if value is None:
            return f'<c r="{reference}"/>'
        escaped = value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        return f'<c r="{reference}" t="inlineStr"><is><t>{escaped}</t></is></c>'

    rows = []
    for row_index, row in enumerate(grid, start=1):
        cells = "".join(
            cell(f"{chr(ord('A') + column)}{row_index}", value)
            for column, value in enumerate(row)
        )
        rows.append(f'<row r="{row_index}">{cells}</row>')
    sheet = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        f'<worksheet xmlns="{_MAIN_NS}"><sheetData>{"".join(rows)}</sheetData></worksheet>'
    )
    workbook = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        f'<workbook xmlns="{_MAIN_NS}" xmlns:r="{_DOC_NS}"><sheets>'
        f'<sheet name="{sheet_name}" sheetId="1" r:id="rId1"/>'
        "</sheets></workbook>"
    )
    rels = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        f'<Relationships xmlns="{_PKG_NS}">'
        '<Relationship Id="rId1" Type="" Target="worksheets/sheet1.xml"/>'
        "</Relationships>"
    )
    with zipfile.ZipFile(path, mode="w") as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", rels)
        archive.writestr("xl/worksheets/sheet1.xml", sheet)


class GSE296875DonorSupplementTests(unittest.TestCase):
    def setUp(self) -> None:
        with TASKSPEC.open("rb") as handle:
            self.spec = tomllib.load(handle)

    def test_taskspec_pins_one_workbook_by_hash_and_size(self) -> None:
        self.assertEqual(len(self.spec["source_supplement_sha256"]), 64)
        self.assertEqual(self.spec["source_supplement_size_bytes"], 5_884_547)
        self.assertEqual(self.spec["source_sheet"], supplement.EXPECTED_SHEET)
        self.assertTrue(
            self.spec["source_supplement"].endswith("-mmc2.xlsx"),
            self.spec["source_supplement"],
        )

    def test_refuses_a_workbook_whose_hash_differs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "mmc2.xlsx"
            _workbook(path, supplement.EXPECTED_SHEET, [["donor_id"], ["1"]])
            with self.assertRaisesRegex(supplement.SupplementAcquisitionError, "SHA-256"):
                supplement.run(workbook=path, taskspec=TASKSPEC, output=root / "out")
            self.assertFalse((root / "out").exists())

    def test_refuses_a_workbook_whose_size_differs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "mmc2.xlsx"
            _workbook(path, supplement.EXPECTED_SHEET, [["donor_id"], ["1"]])
            spec = root / "spec.toml"
            digest = supplement.sha256_file(path)
            text = TASKSPEC.read_text(encoding="utf-8")
            text = text.replace(self.spec["source_supplement_sha256"], digest)
            spec.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(supplement.SupplementAcquisitionError, "byte count"):
                supplement.run(workbook=path, taskspec=spec, output=root / "out")

    def test_refuses_a_workbook_without_the_authoritative_sheet(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "mmc2.xlsx"
            _workbook(path, "S2 other", [["donor_id"], ["1"]])
            spec = self._matching_spec(root, path)
            with self.assertRaisesRegex(supplement.SupplementAcquisitionError, "sheet is absent"):
                supplement.run(workbook=path, taskspec=spec, output=root / "out")

    def test_freezes_a_lossless_sheet_inventory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "mmc2.xlsx"
            _workbook(
                path,
                supplement.EXPECTED_SHEET,
                [
                    ["donor_id", "fibrosis_status", "steatosis_numeric", "BMI"],
                    ["331", "None", "5", "27.4"],
                    ["342", "mild", None, "31.0"],
                    ["346", None, "40", None],
                ],
            )
            spec = self._matching_spec(root, path)
            output = root / "out"
            receipt = supplement.run(workbook=path, taskspec=spec, output=output)
            verify_frozen_tree(output)
            self.assertEqual(receipt["data_rows"], 3)
            schema = json.loads(
                (output / "inventory/s1_donorinfo_schema.json").read_text(encoding="utf-8")
            )
            by_name = {column["name"]: column for column in schema["columns"]}
            self.assertEqual(by_name["fibrosis_status"]["non_blank_rows"], 2)
            self.assertEqual(by_name["fibrosis_status"]["blank_rows"], 1)
            self.assertEqual(by_name["steatosis_numeric"]["non_blank_rows"], 2)
            self.assertEqual(by_name["BMI"]["non_blank_rows"], 2)
            self.assertEqual(
                by_name["fibrosis_status"]["distinct_values"], ["None", "mild"]
            )
            table = (output / "inventory/s1_donorinfo.tsv").read_text(encoding="utf-8")
            self.assertEqual(table.splitlines()[0].split("\t")[0], "donor_id")
            self.assertEqual(table.splitlines()[2], "342\tmild\t\t31.0")

    def test_a_blank_cell_is_never_rendered_as_a_value(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "mmc2.xlsx"
            _workbook(
                path,
                supplement.EXPECTED_SHEET,
                [["donor_id", "steatosis_numeric"], ["331", None]],
            )
            spec = self._matching_spec(root, path)
            output = root / "out"
            supplement.run(workbook=path, taskspec=spec, output=output)
            table = (output / "inventory/s1_donorinfo.tsv").read_text(encoding="utf-8")
            self.assertEqual(table.splitlines()[1], "331\t")
            self.assertNotIn("0\t0", table)

    def test_a_one_cell_prose_caption_is_not_mistaken_for_the_header(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "mmc2.xlsx"
            _workbook(
                path,
                supplement.EXPECTED_SHEET,
                [
                    ["Table S1: Donor characteristics. sample_id is the donor.", None, None],
                    ["sample_id", "age_in_yr", "steatosis_numeric"],
                    ["331", "42", "5"],
                    ["342", "13", "0"],
                ],
            )
            spec = self._matching_spec(root, path)
            output = root / "out"
            receipt = supplement.run(workbook=path, taskspec=spec, output=output)
            self.assertEqual(receipt["header"], ["sample_id", "age_in_yr", "steatosis_numeric"])
            self.assertEqual(receipt["header_row_index"], 1)
            self.assertEqual(receipt["data_rows"], 2)
            table = (output / "inventory/s1_donorinfo.tsv").read_text(encoding="utf-8")
            self.assertEqual(
                table.splitlines(),
                ["sample_id\tage_in_yr\tsteatosis_numeric", "331\t42\t5", "342\t13\t0"],
            )

    def test_the_na_sentinel_is_counted_as_missing_not_observed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "mmc2.xlsx"
            _workbook(
                path,
                supplement.EXPECTED_SHEET,
                [
                    ["sample_id", "steatosis_numeric"],
                    ["331", "5"],
                    ["342", "NA"],
                    ["346", None],
                ],
            )
            spec = self._matching_spec(root, path)
            output = root / "out"
            supplement.run(workbook=path, taskspec=spec, output=output)
            schema = json.loads(
                (output / "inventory/s1_donorinfo_schema.json").read_text(encoding="utf-8")
            )
            column = {item["name"]: item for item in schema["columns"]}["steatosis_numeric"]
            self.assertEqual(column["non_blank_rows"], 2)
            self.assertEqual(column["sentinel_missing_rows"], 1)
            self.assertEqual(column["blank_rows"], 1)
            self.assertEqual(column["observed_rows"], 1)
            self.assertEqual(schema["missing_sentinels"], ["NA"])

    def _matching_spec(self, root: Path, workbook: Path) -> Path:
        spec = root / "spec.toml"
        text = TASKSPEC.read_text(encoding="utf-8")
        text = text.replace(
            self.spec["source_supplement_sha256"], supplement.sha256_file(workbook)
        )
        text = text.replace(
            f'source_supplement_size_bytes = {self.spec["source_supplement_size_bytes"]}',
            f"source_supplement_size_bytes = {workbook.stat().st_size}",
        )
        spec.write_text(text, encoding="utf-8")
        return spec


if __name__ == "__main__":
    unittest.main()
