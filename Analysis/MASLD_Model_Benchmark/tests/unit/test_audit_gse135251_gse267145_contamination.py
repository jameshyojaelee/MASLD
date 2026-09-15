from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.audit_gse135251_gse267145_contamination import (
    ContaminationAuditError,
    describe,
    per_sample_rank,
    read_rowname_offset_tsv,
    scan_namespaces,
    upper_triangle,
)


class RownameOffsetTests(unittest.TestCase):
    """The N+1 defect shifts every column; it must fail closed, not silently."""

    def test_valid_n_plus_one_file_parses(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "meta.tsv"
            path.write_text(
                "sample_id\tRun\tnas_score\n"
                "SRR1\tSRR1\tSRR1\t4\n"
                "SRR2\tSRR2\tSRR2\t7\n",
                encoding="utf-8",
            )
            header, rows = read_rowname_offset_tsv(path)
            self.assertEqual(header, ["sample_id", "Run", "nas_score"])
            # header[i] maps to row[i+1]; row[0] is the unnamed rowname.
            self.assertEqual(rows[0][header.index("nas_score") + 1], "4")
            self.assertEqual(rows[1][header.index("nas_score") + 1], "7")

    def test_a_well_formed_file_without_rownames_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "meta.tsv"
            path.write_text("sample_id\tRun\tnas_score\nSRR1\tSRR1\t4\n", encoding="utf-8")
            with self.assertRaises(ContaminationAuditError):
                read_rowname_offset_tsv(path)

    def test_ragged_row_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "meta.tsv"
            path.write_text(
                "a\tb\tc\nr1\t1\t2\t3\nr2\t1\t2\n", encoding="utf-8"
            )
            with self.assertRaises(ContaminationAuditError):
                read_rowname_offset_tsv(path)


class NamespaceScanTests(unittest.TestCase):
    def test_every_namespace_is_found(self) -> None:
        text = "GSM3998167 SRX6635477 SRR9882956 SAMN12426403 PRJNA558102 SRP217231"
        found = scan_namespaces(text)
        self.assertEqual(found["geo_sample"], {"GSM3998167"})
        self.assertEqual(found["sra_experiment"], {"SRX6635477"})
        self.assertEqual(found["sra_run"], {"SRR9882956"})
        self.assertEqual(found["biosample"], {"SAMN12426403"})
        self.assertEqual(found["bioproject"], {"PRJNA558102"})
        self.assertEqual(found["sra_study"], {"SRP217231"})

    def test_accession_prefixes_do_not_bleed(self) -> None:
        found = scan_namespaces("GSE135251 GSM1")
        self.assertEqual(found["geo_sample"], {"GSM1"})


class RankSimilarityTests(unittest.TestCase):
    def test_rank_correlation_is_scale_free(self) -> None:
        """A duplicate sequenced to a different depth must still look identical."""

        rng = np.random.default_rng(5)
        profile = rng.lognormal(size=(1, 500))
        deep = profile * 17.0
        ranked = per_sample_rank(np.vstack([profile, deep]))
        self.assertAlmostEqual(float(ranked[0] @ ranked[1]), 1.0, places=12)

    def test_independent_profiles_are_near_zero(self) -> None:
        rng = np.random.default_rng(9)
        ranked = per_sample_rank(rng.lognormal(size=(2, 5000)))
        self.assertLess(abs(float(ranked[0] @ ranked[1])), 0.1)

    def test_a_planted_duplicate_is_the_reciprocal_best_hit(self) -> None:
        rng = np.random.default_rng(17)
        source = rng.lognormal(size=(12, 800))
        target = rng.lognormal(size=(9, 800))
        # Plant source row 5 into target row 3 at a different depth.
        target[3] = source[5] * 3.5
        cross = per_sample_rank(source) @ per_sample_rank(target).T
        self.assertEqual(int(np.argmax(cross, axis=1)[5]), 3)
        self.assertEqual(int(np.argmax(cross, axis=0)[3]), 5)
        self.assertGreater(float(cross[5, 3]), 0.99)

    def test_flat_sample_is_rejected(self) -> None:
        with self.assertRaises(ContaminationAuditError):
            per_sample_rank(np.ones((2, 10)))


class SummaryTests(unittest.TestCase):
    def test_upper_triangle_excludes_the_diagonal(self) -> None:
        matrix = np.asarray([[1.0, 0.2, 0.3], [0.2, 1.0, 0.4], [0.3, 0.4, 1.0]])
        self.assertEqual(sorted(upper_triangle(matrix).tolist()), [0.2, 0.3, 0.4])

    def test_describe_reports_the_pair_count(self) -> None:
        stats = describe(np.arange(100, dtype=np.float64))
        self.assertEqual(stats["n_pairs"], 100)
        self.assertAlmostEqual(stats["max"], 99.0)


if __name__ == "__main__":
    unittest.main()
