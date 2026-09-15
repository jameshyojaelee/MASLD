from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.overlay_lineage_variant_provenance import (
    BACKUP,
    CANONICAL,
    classify_quantification,
    counts_fingerprint,
)


class GlobMechanismTests(unittest.TestCase):
    """The defect this overlay records: a filter that reads like a listing."""

    def test_the_canonical_glob_does_not_match_the_backup(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            root = Path(value) / "GSE1"
            root.mkdir()
            (root / f"GSE1{CANONICAL}").write_text("x\n", encoding="utf-8")
            (root / f"GSE1{BACKUP}").write_text("x\n", encoding="utf-8")
            matched = sorted(p.name for p in root.glob(f"*{CANONICAL}"))
            self.assertEqual(matched, [f"GSE1{CANONICAL}"])
            self.assertNotIn(f"GSE1{BACKUP}", matched)
            self.assertEqual(len(list(root.iterdir())), 2, "both files exist on disk")

    def test_the_suffixes_are_distinct_not_nested(self) -> None:
        self.assertFalse(BACKUP.endswith(CANONICAL))
        self.assertTrue(BACKUP.startswith(CANONICAL[:-4]))


class QuantificationClassificationTests(unittest.TestCase):
    def test_gene_universes_separate_the_pipelines(self) -> None:
        self.assertEqual(classify_quantification(77078, "TSPAN6"), "kallisto_human")
        self.assertEqual(classify_quantification(37606, "A1BG"), "star_human")
        self.assertEqual(classify_quantification(37615, "A1BG"), "star_human")
        self.assertEqual(classify_quantification(33423, "0610005C13Rik"), "mouse")

    def test_an_unrecognised_universe_is_named_not_guessed(self) -> None:
        out = classify_quantification(12345, "FOO")
        self.assertTrue(out.startswith("unrecognised_"))
        self.assertIn("12345", out)


class CountsFingerprintTests(unittest.TestCase):
    def test_rows_and_columns_are_counted_without_loading_the_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            p = Path(value) / "counts.tsv"
            p.write_text(
                "SYMBOL\tSRR1\tSRR2\tSRR3\n"
                "TSPAN6\t1\t2\t3\n"
                "A1BG\t4\t5\t6\n"
                "TP53\t7\t8\t9\n",
                encoding="utf-8",
            )
            fp = counts_fingerprint(p)
            self.assertEqual(fp["gene_rows"], 3)
            self.assertEqual(fp["sample_columns"], 3)
            self.assertEqual(fp["first_row_key"], "TSPAN6")

    def test_trailing_blank_lines_do_not_inflate_the_gene_count(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            p = Path(value) / "counts.tsv"
            p.write_text("SYMBOL\tS1\nTSPAN6\t1\nA1BG\t2\n\n\n", encoding="utf-8")
            self.assertEqual(counts_fingerprint(p)["gene_rows"], 2)


if __name__ == "__main__":
    unittest.main()
