"""Fixture tests for step-internal helpers (peak index, gene-version handling, repeated donors)."""

from __future__ import annotations

import importlib.util
import os
import pathlib
import sys
import tempfile
import unittest

import anndata
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
os.environ.setdefault("AGA_OUT_ROOT", tempfile.mkdtemp())


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), HERE / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestPeakIndex(unittest.TestCase):
    def test_variant_inside_peak_hits_once_and_outside_misses(self):
        m = load("10_package_B_chromatin.py")
        idx = m.PeakIndex()
        idx.add("hepatocyte", "chr1", 100, 600, "p1")      # 0-based half-open: covers 1-based 101..600
        idx.add("hepatocyte", "chr1", 1000, 1500, "p2")
        idx.build()
        self.assertEqual(idx.hits("hepatocyte", "chr1", 101), ["p1"])
        self.assertEqual(idx.hits("hepatocyte", "chr1", 600), ["p1"])
        self.assertEqual(idx.hits("hepatocyte", "chr1", 100), [])       # window/peak mismatch: position before the peak start
        self.assertEqual(idx.hits("hepatocyte", "chr1", 601), [])
        self.assertEqual(idx.hits("stellate", "chr1", 300), [])         # missing lineage returns nothing, no fallback

    def test_gc_content(self):
        m = load("10_package_B_chromatin.py")

        class F:
            def fetch(self, chrom, s, e):
                return "GGCCAATT"[s:e]

        self.assertAlmostEqual(m.gc_content(F(), "chr1:0-8"), 0.5)


class TestObsFrame(unittest.TestCase):
    def test_gene_version_is_stripped_and_missing_columns_filled(self):
        m = load("07_build_prediction_records.py")
        obs = pd.DataFrame({"variant": ["chr1:5:A>G"], "gene_id": ["ENSG00000134243.12"], "gene_name": ["SORT1"], "strand": ["+"]}, index=["0"])
        a = anndata.AnnData(X=np.zeros((1, 1), dtype=np.float32), obs=obs, var=pd.DataFrame({"name": ["t"]}, index=["0"]))
        f = m.obs_frame(a)
        self.assertEqual(f.at[0, "gene_id"], "ENSG00000134243")
        self.assertEqual(f.at[0, "variant_uid"], "chr1:5:A:G")
        self.assertEqual(f.at[0, "junction_Start"], "")

    def test_panel_group_prefers_ontology_and_falls_back_to_hepg2_name(self):
        m = load("07_build_prediction_records.py")
        self.assertEqual(m.panel_group("UBERON:0002107", "liver"), "primary_liver")
        self.assertEqual(m.panel_group("EFO:9999999", "HepG2 cells"), "HepG2")
        self.assertIsNone(m.panel_group("UBERON:0002048", "lung"))


class TestRepeatedDonors(unittest.TestCase):
    def test_measured_context_uses_recorded_donor_counts_not_row_counts(self):
        # da_lineage_summary carries donor denominators; duplicating rows in a derived table must not change them
        rows = [{"lineage": "stellate", "n_normal_gse244832": "5", "n_mash_gse244832": "9"}] * 3
        donors = {(r["lineage"], r["n_normal_gse244832"], r["n_mash_gse244832"]) for r in rows}
        self.assertEqual(donors, {("stellate", "5", "9")})


if __name__ == "__main__":
    unittest.main()


class TestBenchmarkUids(unittest.TestCase):
    def test_validated_uid_checks_reference_and_classifies_indels(self):
        m = load("31_query_benchmark_sets.py")

        class Fasta:
            references = ["chr1"]

            def fetch(self, chrom, s, e):
                return "GATTACA"[s:e]
        f = Fasta(); stats = {"snv": 0, "indel_not_served": 0, "fasta_mismatch": 0, "nonstandard_chrom": 0}
        self.assertEqual(m.validated_uid(f, "1", 2, "A", "G", stats), "chr1:2:A:G")          # chr prefix added, ref A at 1-based 2
        self.assertIsNone(m.validated_uid(f, "chr1", 3, "TT", "T", stats))                  # indel: reference checked, then NOT served (probed 2026-09-09)
        self.assertIsNone(m.validated_uid(f, "chr1", 3, "TA", "T", stats))                  # second base mismatch -> None, never repaired
        self.assertIsNone(m.validated_uid(f, "chr1", 2, "C", "G", stats))                   # SNV mismatch
        self.assertIsNone(m.validated_uid(f, "chrUn_1", 2, "A", "G", stats))
        self.assertEqual(stats, {"snv": 1, "indel_not_served": 1, "fasta_mismatch": 2, "nonstandard_chrom": 1})


class TestObsFrameEmpty(unittest.TestCase):
    def test_scorer_with_no_rows_yields_empty_frame_with_expected_columns(self):
        m = load("07_build_prediction_records.py")
        empty = anndata.AnnData(X=np.zeros((0, 3), dtype=np.float32), var=pd.DataFrame(index=["t1", "t2", "t3"]))
        obs = m.obs_frame(empty)                       # a chunk can carry zero variants for a scorer (e.g. no splice junction nearby)
        self.assertEqual(len(obs), 0)
        self.assertEqual(list(obs.columns), ["variant_uid", "gene_id", "gene_name", "strand", "junction_Start", "junction_End"])


class TestWriteOnceDir(unittest.TestCase):
    """Step 14's guard must refuse a directory that HOLDS results, not an empty one made by setup."""

    def test_empty_existing_dir_is_accepted_and_nonempty_is_refused(self):
        import lib_atlas as la
        m = load("14_validate.py")
        with tempfile.TemporaryDirectory() as d:
            root = pathlib.Path(d)
            fresh = root / "fresh"
            m.refuse_if_written(fresh)                 # does not exist -> created
            self.assertTrue(fresh.is_dir())
            m.refuse_if_written(fresh)                 # exists but empty -> accepted
            (fresh / "validation.json").write_text("{}")
            with self.assertRaises(la.ContractError):  # holds a result -> refused
                m.refuse_if_written(fresh)
