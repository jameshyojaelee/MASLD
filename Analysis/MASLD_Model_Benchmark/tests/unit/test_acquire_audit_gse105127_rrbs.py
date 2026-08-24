from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts.acquire_audit_gse105127_rrbs import (
    https_url,
    parse_awk_output,
    read_rrbs_rows,
)


class GSE105127RRBSAcquisitionTests(unittest.TestCase):
    def test_ncbi_ftp_url_is_promoted_to_https(self) -> None:
        observed = https_url("ftp://ftp.ncbi.nlm.nih.gov/a/file.CG.bed.gz")
        self.assertEqual(observed, "https://ftp.ncbi.nlm.nih.gov/a/file.CG.bed.gz")

    def test_awk_output_preserves_large_integer_counts(self) -> None:
        observed = parse_awk_output(
            "rows\t14000000\nweighted_methylated_sum\t12.500000\nerror_total\t0\n"
        )
        self.assertEqual(observed["rows"], 14_000_000)
        self.assertEqual(observed["weighted_methylated_sum"], 12.5)

    def test_join_reader_rejects_fewer_than_57_rrbs_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample_join.tsv").write_text(
                "assay\tsample_accession\tpairing\nRRBS\tGSM1\tadjacent_section\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(Exception, "57 unique RRBS"):
                read_rrbs_rows(root)


if __name__ == "__main__":
    unittest.main()
