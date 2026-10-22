"""Coordinate invariants; run on a compute node with the existing figure env.

Set FIG2_PROJECT_ROOT to enable the fixed-source scientific regression.
"""
import os
from pathlib import Path
import unittest

import pandas as pd

from fig2_locus_coordinates import recover_locus_positions


class CoordinateRecoveryTests(unittest.TestCase):
    def recover(self, locus, positions, chromosome="2", **kwargs):
        frame = pd.DataFrame({"credible_set_uid": ["cs"], "locus": [locus]})
        members = pd.DataFrame({"credible_set_uid": ["cs"] * len(positions),
                                "chromosome": [chromosome] * len(positions),
                                "position": positions})
        return recover_locus_positions(frame, members, **kwargs)

    def test_gckr_and_float_encoding(self):
        for locus in ["2.2773094", 2.2773094]:
            with self.subTest(locus=locus):
                row = self.recover(locus, [27730940, 27731000]).iloc[0]
                self.assertEqual(row["pos"], 27730940)
                self.assertEqual(row["encoded_pos"], 2773094)
                self.assertEqual(row["trailing_zeros_restored"], 1)

    def test_unchanged_and_multiple_zero_restoration(self):
        for locus, position, zeros, chrom in [
                ("22.44324727", 44324727, 0, "22"),
                ("2.277309", 27730900, 2, "2"),
                ("2.1000", 1000, 0, "2")]:
            with self.subTest(locus=locus):
                row = self.recover(locus, [position], chromosome=chrom,
                                   max_distance_bp=0).iloc[0]
                self.assertEqual(row["pos"], position)
                self.assertEqual(row["trailing_zeros_restored"], zeros)

    def test_sentinel_need_not_be_member_or_mean(self):
        row = self.recover("2.2773094", [27730920, 27730990],
                           max_distance_bp=100).iloc[0]
        self.assertEqual(row["pos"], 27730940)
        self.assertNotEqual(row["pos"], (27730920 + 27730990) / 2)

    def test_all_member_endpoints_bound_recovery(self):
        self.assertEqual(self.recover("2.2773094", [27730840, 27731040],
                                      max_distance_bp=100).iloc[0]["pos"], 27730940)
        with self.assertRaises(ValueError):
            self.recover("2.2773094", [27730840, 27731041], max_distance_bp=100)

    def test_ambiguity_and_unrecoverable_coordinates_fail(self):
        for locus, positions in [("2.1", [1, 10]), ("2.01", [1, 10]),
                                 ("2.2773094", [80000000])]:
            with self.subTest(locus=locus), self.assertRaises(ValueError):
                self.recover(locus, positions)

    def test_native_coordinates_must_be_positive_integers(self):
        for position in [0, -1, 27730940.5, None, "bad"]:
            with self.subTest(position=position), self.assertRaises(ValueError):
                self.recover("2.2773094", [position])

    def test_invalid_locus_fails(self):
        for locus in [None, "", "2", "2.bad", "2.0", "0.123", "2.12.3"]:
            with self.subTest(locus=locus), self.assertRaises(ValueError):
                self.recover(locus, [27730940])

    def test_chromosome_disagreement_and_missing_uid_fail(self):
        f = pd.DataFrame({"credible_set_uid": ["a"], "locus": ["2.2773094"]})
        for m in [pd.DataFrame({"credible_set_uid": ["a", "a"],
                               "chromosome": ["2", "3"], "position": [27730940]*2}),
                  pd.DataFrame({"credible_set_uid": ["a"], "chromosome": ["3"],
                                "position": [27730940]}),
                  pd.DataFrame({"credible_set_uid": ["b"], "chromosome": ["2"],
                                "position": [27730940]})]:
            with self.assertRaises(ValueError):
                recover_locus_positions(f, m)

    def test_order_index_and_sources_preserved(self):
        f = pd.DataFrame({"credible_set_uid": ["b", "a"],
                          "locus": ["22.44324727", "2.2773094"],
                          "mass": [0.2, 0.8]}, index=[17, 3])
        m = pd.DataFrame({"credible_set_uid": ["a", "b", "a"],
                          "chromosome": ["2", "22", "2"],
                          "position": [27731000, 44324727, 27730940]})
        before_f, before_m = f.copy(deep=True), m.copy(deep=True)
        result = recover_locus_positions(f, m)
        pd.testing.assert_frame_equal(f, before_f)
        pd.testing.assert_frame_equal(m, before_m)
        pd.testing.assert_frame_equal(result[f.columns], f)
        self.assertEqual(result["pos"].tolist(), [44324727, 27730940])
        self.assertEqual(result["chrom"].tolist(), ["22", "2"])
        self.assertEqual(result["member_min_position"].tolist(), [44324727, 27730940])
        self.assertEqual(result["member_max_position"].tolist(), [44324727, 27731000])

    def test_optional_member_locus_cannot_disagree(self):
        f = pd.DataFrame({"credible_set_uid": ["a"], "locus": ["2.2773094"]})
        m = pd.DataFrame({"credible_set_uid": ["a"], "chromosome": ["2"],
                          "position": [27730940], "locus": ["2.8888888"]})
        with self.assertRaises(ValueError):
            recover_locus_positions(f, m)

    def test_missing_columns_and_empty_inputs_fail(self):
        f = pd.DataFrame({"credible_set_uid": ["a"], "locus": ["2.2773094"]})
        m = pd.DataFrame({"credible_set_uid": ["a"], "chromosome": ["2"],
                          "position": [27730940]})
        cases = [(f.iloc[:0], m), (f, m.iloc[:0])]
        cases += [(f.drop(columns=c), m) for c in f.columns]
        cases += [(f, m.drop(columns=c)) for c in m.columns]
        for a, b in cases:
            with self.subTest(frame_columns=list(a), member_columns=list(b)):
                with self.assertRaises((ValueError, KeyError)):
                    recover_locus_positions(a, b)


@unittest.skipUnless(os.environ.get("FIG2_PROJECT_ROOT"), "fixed-source regression disabled")
class FixedSourceRegression(unittest.TestCase):
    def test_gated_physical_cluster_coding_contrast(self):
        from fig2_multiplicity_correction import apply_gate, ARCH_COLS, SCOPES
        base = Path(os.environ["FIG2_PROJECT_ROOT"]) / "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321"
        a = pd.read_csv(base / "credible_set_pip_architecture.tsv", sep="\t", dtype={"locus": str})
        m = pd.read_csv(base / "credible_set_members_annotated.tsv", sep="\t", dtype={"locus": str}, low_memory=False)
        a = a[a["trait_scope"].isin(SCOPES)].copy()
        m = m[m["trait_scope"].isin(SCOPES)].copy()
        a, m, _ = apply_gate(a, m)
        r = recover_locus_positions(a, m)
        self.assertEqual(len(r), 6415)
        self.assertEqual(int((r["trailing_zeros_restored"] > 0).sum()), 715)
        # Independent sorted adjacency clustering, with no producer cluster call.
        groups = {}
        last_chrom, last_pos, number = None, None, -1
        for row in sorted(r.to_dict("records"), key=lambda x: (x["chrom"], x["pos"])):
            if row["chrom"] != last_chrom or row["pos"] - last_pos > 1000000:
                number += 1
            groups.setdefault((number, row["trait_scope"]), []).append(row)
            last_chrom, last_pos = row["chrom"], row["pos"]
        direct = {g for g, s in groups if s == SCOPES[0]}
        enzyme = {g for g, s in groups if s == SCOPES[1]}
        self.assertEqual((len(direct | enzyme), len(direct), len(enzyme), len(direct & enzyme)),
                         (233, 19, 228, 14))
        fractions = []
        for scope in SCOPES:
            numerator = denominator = 0.0
            for (_, s), rows in groups.items():
                if s == scope:
                    numerator += sum(x["protein_altering_pip_mass"] + x["canonical_splice_pip_mass"] for x in rows) / len(rows)
                    denominator += sum(sum(x[c] for c in ARCH_COLS) for x in rows) / len(rows)
            fractions.append(numerator / denominator)
        self.assertAlmostEqual(fractions[0] / fractions[1], 6.026671, places=5)
        self.assertAlmostEqual(fractions[0], 0.2334, delta=0.0001)
        self.assertAlmostEqual(fractions[1], 0.0387, delta=0.0001)
        print("Fixed-source regression: 6415 sets, 715 restored; clusters 233/19/228/14; fractions", fractions)


if __name__ == "__main__":
    unittest.main()
