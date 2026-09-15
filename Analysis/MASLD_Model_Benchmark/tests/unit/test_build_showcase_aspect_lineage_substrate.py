from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts.build_showcase_aspect_lineage_substrate import (
    REF_HUMAN_LINEAGES,
    SAF_TRIPLE,
    SubstrateError,
    build_shared_gene_axis,
    parse_series_matrix_characteristics,
    read_bayesprism_proportions,
    read_flat_tsv,
    read_sra_run_info,
)


class FlatTableReaderTests(unittest.TestCase):
    """A flat reader must refuse an R rowname-offset table, not shift columns."""

    def test_a_flat_table_reads(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "flat.tsv"
            path.write_text(
                "sample_id\tdisease\tFibrosis_stage\n"
                "SRR1\tNAFLD\t2\n",
                encoding="utf-8",
            )
            header, rows = read_flat_tsv(path)
            self.assertEqual(header, ["sample_id", "disease", "Fibrosis_stage"])
            self.assertEqual(rows[0]["Fibrosis_stage"], "2")

    def test_a_rowname_offset_table_fails_closed(self) -> None:
        """The naive read returns the neighbouring column and never raises."""

        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "offset.tsv"
            path.write_text(
                "sample_id\tdisease\tFibrosis_stage\n"
                "SRR1\tSRR1\tNAFLD\t2\n",
                encoding="utf-8",
            )
            with self.assertRaises(SubstrateError):
                read_flat_tsv(path)


class BayesPrismReaderTests(unittest.TestCase):
    def _write(self, path: Path, *, with_key: bool = True) -> None:
        lines = ["\t".join(REF_HUMAN_LINEAGES)]
        row = [f"0.{index:04d}" for index in range(len(REF_HUMAN_LINEAGES))]
        total = sum(float(v) for v in row)
        row[0] = f"{float(row[0]) + 1.0 - total:.10f}"
        lines.append("\t".join((["SRR1"] if with_key else []) + row))
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def test_the_unnamed_key_column_is_not_folded_into_a_lineage(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "prop.tsv"
            self._write(path)
            lineages, keys, matrix = read_bayesprism_proportions(path)
            self.assertEqual(tuple(lineages), REF_HUMAN_LINEAGES)
            self.assertEqual(keys, ["SRR1"])
            self.assertEqual(matrix.shape, (1, 16))
            self.assertAlmostEqual(float(matrix.sum()), 1.0, places=9)

    def test_sixteen_lineages_not_fifteen(self) -> None:
        self.assertEqual(len(REF_HUMAN_LINEAGES), 16)

    def test_a_table_without_the_key_column_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "prop.tsv"
            self._write(path, with_key=False)
            with self.assertRaises(SubstrateError):
                read_bayesprism_proportions(path)


class SeriesMatrixParsingTests(unittest.TestCase):
    def _write(self, path: Path) -> None:
        # Row three is positionally ragged on purpose: it carries liver lobe for
        # the first sample and gender for the second.  Keying the row by its
        # first cell would file a gender under liver lobe.
        lines = [
            '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"',
            '!Sample_characteristics_ch1\t"patient id: P1"\t"patient id: P2"'
            '\t"patient id: P1"',
            '!Sample_characteristics_ch1\t"liver lobe: Right"\t"gender: F"\t""',
            '!Sample_characteristics_ch1\t"saf score: S2A3F3"'
            '\t"saf score: end stage"\t"saf score: S2A3F3"',
        ]
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")

    def test_every_cell_is_split_on_its_own_key(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "series.txt.gz"
            self._write(path)
            accessions, fields = parse_series_matrix_characteristics(path)
            self.assertEqual(accessions, ["GSM1", "GSM2", "GSM3"])
            self.assertEqual(fields["GSM1"]["liver_lobe"], "Right")
            self.assertEqual(fields["GSM2"]["gender"], "F")
            self.assertNotIn("gender", fields["GSM1"])
            self.assertNotIn("liver_lobe", fields["GSM2"])

    def test_the_saf_field_is_present_even_when_ungraded(self) -> None:
        """Ungraded is a failed pattern match, never an absent field."""

        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "series.txt.gz"
            self._write(path)
            _, fields = parse_series_matrix_characteristics(path)
            self.assertEqual(fields["GSM2"]["saf_score"], "end stage")
            self.assertIsNone(SAF_TRIPLE.fullmatch(fields["GSM2"]["saf_score"]))
            self.assertIsNotNone(SAF_TRIPLE.fullmatch(fields["GSM1"]["saf_score"]))

    def test_two_gsms_of_one_donor_do_not_become_two_donors(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "series.txt.gz"
            self._write(path)
            accessions, fields = parse_series_matrix_characteristics(path)
            self.assertEqual(len(accessions), 3)
            self.assertEqual(len({fields[a]["patient_id"] for a in accessions}), 2)

    def test_conflicting_values_for_one_key_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "series.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(
                    '!Sample_geo_accession\t"GSM1"\n'
                    '!Sample_characteristics_ch1\t"patient id: P1"\n'
                    '!Sample_characteristics_ch1\t"patient id: P2"\n'
                )
            with self.assertRaises(SubstrateError):
                parse_series_matrix_characteristics(path)


class SraRunInfoTests(unittest.TestCase):
    def test_columns_are_located_by_content_on_a_headerless_file(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "SraRunInfo.csv"
            path.write_text(
                'SRR1,x,SRX1,GSM1,"A CENTRE\n NAME",GSM1\n'
                'SRR2,y,SRX2,GSM1,"A CENTRE\n NAME",GSM1\n',
                encoding="utf-8",
            )
            rows = read_sra_run_info(path)
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["run_accession"], "SRR1")
            self.assertEqual(rows[0]["experiment_accession"], "SRX1")
            self.assertEqual(rows[0]["sample_accession"], "GSM1")
            self.assertEqual(
                len({row["sample_accession"] for row in rows}),
                1,
                "two runs of one sample must not become two samples",
            )

    def test_duplicate_runs_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "SraRunInfo.csv"
            path.write_text("SRR1,SRX1,GSM1\nSRR1,SRX2,GSM2\n", encoding="utf-8")
            with self.assertRaises(SubstrateError):
                read_sra_run_info(path)


class SharedGeneAxisTests(unittest.TestCase):
    def test_intersection_indices_point_at_the_same_gene(self) -> None:
        a = ["ENSG00000000003", "ENSG00000000005", "ENSG00000000419"]
        b = ["ENSG00000000419", "ENSG00000000003", "ENSG00000001036"]
        shared, a_cols, b_cols = build_shared_gene_axis(a, b)
        self.assertEqual(shared, ["ENSG00000000003", "ENSG00000000419"])
        for index, gene in enumerate(shared):
            self.assertEqual(a[a_cols[index]], gene)
            self.assertEqual(b[b_cols[index]], gene)

    def test_restriction_reorders_both_matrices_consistently(self) -> None:
        a = ["ENSG00000000003", "ENSG00000000005", "ENSG00000000419"]
        b = ["ENSG00000000419", "ENSG00000000003", "ENSG00000001036"]
        shared, a_cols, b_cols = build_shared_gene_axis(a, b)
        left = np.arange(6, dtype=np.float64).reshape(2, 3)
        right = np.arange(100, 106, dtype=np.float64).reshape(2, 3)
        self.assertEqual(list(left[:, a_cols][0]), [0.0, 2.0])
        self.assertEqual(list(right[:, b_cols][0]), [101.0, 100.0])
        self.assertEqual(len(shared), 2)

    def test_versioned_ids_fail_closed(self) -> None:
        with self.assertRaises(SubstrateError):
            build_shared_gene_axis(["ENSG00000000003.15"], ["ENSG00000000003"])

    def test_duplicate_ids_fail_closed(self) -> None:
        with self.assertRaises(SubstrateError):
            build_shared_gene_axis(
                ["ENSG00000000003", "ENSG00000000003"], ["ENSG00000000003"]
            )

    def test_a_symbol_axis_is_not_silently_accepted(self) -> None:
        """The column labelled SYMBOL is 46.3% ENSG; names are not evidence."""

        with self.assertRaises(SubstrateError):
            build_shared_gene_axis(["TP53", "ENSG00000000003"], ["ENSG00000000003"])


if __name__ == "__main__":
    unittest.main()
