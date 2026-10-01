"""Fixture tests for the Figure 2 coding-mass gate and physical-cluster bootstrap
(`fig2_multiplicity_correction.py`, review items 3 and 4)."""

from __future__ import annotations

import pathlib
import sys
import unittest

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import fig2_multiplicity_correction as m  # noqa: E402

DIRECT, ENZYME = m.SCOPES
N_BOOT = 200


def tab(rows):
    """(trait_scope, cluster_id, num, den) rows as a cluster_table output."""
    return pd.DataFrame(rows, columns=["trait_scope", "cluster_id", "num", "den"])


class TestPhysicalClusterBootstrap(unittest.TestCase):
    def test_shared_rows_move_together(self):
        # Every cluster holds one direct and one identical enzyme row. If the
        # pair moves together, each replicate has direct fraction == enzyme fraction.
        rows = []
        for i, (num, den) in enumerate([(0.9, 1.0), (0.1, 1.0), (0.0, 1.0),
                                        (0.5, 2.0), (0.3, 1.0), (0.7, 1.5)]):
            rows += [(DIRECT, f"c{i}", num, den), (ENZYME, f"c{i}", num, den)]
        wide = m.physical_wide(tab(rows))
        self.assertEqual(len(wide), 6)
        da, db = m.cluster_bootstrap(wide, np.random.default_rng(1), N_BOOT)
        np.testing.assert_array_equal(da, db)

        res = m.analyse(tab(rows), np.random.default_rng(1), n_boot=N_BOOT)
        self.assertEqual(res["n_shared_clusters"], 6)
        self.assertEqual(res["diff_boot_lo"], 0.0)
        self.assertEqual(res["diff_boot_hi"], 0.0)

        # The check can fail: resampling each class's rows independently breaks the pairs.
        g = np.random.default_rng(1)
        num, den = wide["num_a"].to_numpy(), wide["den_a"].to_numpy()
        ia, ib = g.integers(0, 6, (N_BOOT, 6)), g.integers(0, 6, (N_BOOT, 6))
        indep = num[ia].sum(1) / den[ia].sum(1) - num[ib].sum(1) / den[ib].sum(1)
        self.assertTrue((indep != 0).any())

    def test_exclusive_class_contributes_nothing_to_the_other(self):
        rows = [(DIRECT, "s", 0.5, 1.0), (ENZYME, "s", 0.2, 1.0),
                (DIRECT, "d1", 1.0, 1.0), (DIRECT, "d2", 0.0, 1.0),
                (ENZYME, "e1", 0.0, 1.0), (ENZYME, "e2", 0.1, 1.0)]
        res = m.analyse(tab(rows), np.random.default_rng(2), n_boot=N_BOOT)
        self.assertEqual((res["n_clusters_direct"], res["n_clusters_enzyme"],
                          res["n_physical_clusters"], res["n_shared_clusters"]), (3, 3, 5, 1))
        self.assertAlmostEqual(res["fraction_direct"], 1.5 / 3)
        self.assertAlmostEqual(res["fraction_enzyme"], 0.3 / 3)

    def test_no_permutation_p_is_reported(self):
        rows = [(s, f"c{i}", 0.1 * i, 1.0) for i in range(4) for s in (DIRECT, ENZYME)]
        res = m.analyse(tab(rows), np.random.default_rng(3), n_boot=N_BOOT)
        self.assertTrue(np.isnan(res["p_permutation"]))
        self.assertTrue(res["p_permutation_status"].startswith("not_computed"))

    def test_bootstrap_p_is_a_floored_nominal_diagnostic_with_no_call(self):
        # every draw has direct > enzyme, so the empty tail gives the 2/(n+1) floor
        rows = [(DIRECT, f"c{i}", 0.9, 1.0) for i in range(4)] + \
               [(ENZYME, f"c{i}", 0.1, 1.0) for i in range(4)]
        res = m.analyse(tab(rows), np.random.default_rng(4), n_boot=N_BOOT)
        self.assertNotIn("p_bootstrap_diff", res)
        self.assertAlmostEqual(res["nominal_p_bootstrap_diff"], 2 / (N_BOOT + 1))
        self.assertFalse([k for k in res if k.startswith("survives")])

    def test_pooled_positions_give_one_id_across_classes(self):
        frame = pd.DataFrame({"trait_scope": [DIRECT, ENZYME, ENZYME],
                              "locus": ["1.100", "1.500000", "1.5000000"]})
        ids = m.assign_clusters(m.parse_locus(frame)).set_index("locus")["cluster_id"]
        self.assertEqual(ids["1.100"], ids["1.500000"])
        self.assertNotEqual(ids["1.100"], ids["1.5000000"])

    def test_fit_seed_does_not_depend_on_run_order(self):
        first = m.fit_rng("window_size|1Mb").integers(0, 10**9, 5)
        m.fit_rng("reference_1Mb|none").integers(0, 10**9, 5)
        again = m.fit_rng("window_size|1Mb").integers(0, 10**9, 5)
        np.testing.assert_array_equal(first, again)


def fixture(gate, masses, member_pips):
    """Architecture and member tables for credible sets cs0, cs1, ..."""
    arch = pd.DataFrame({
        "credible_set_uid": [f"cs{i}" for i in range(len(gate))],
        "trait_scope": [DIRECT] + [ENZYME] * (len(gate) - 1),
        "study": "s", "pip_sum_gate_passed": gate,
    })
    for j, col in enumerate(m.ARCH_COLS):
        arch[col] = [row[j] for row in masses]
    members = pd.DataFrame(
        [(f"cs{i}", arch.trait_scope[i], "s", gate[i], p)
         for i, pips in enumerate(member_pips) for p in pips],
        columns=["credible_set_uid", "trait_scope", "study",
                 "pip_sum_gate_passed", "normalized_susie_pip"])
    return arch, members


class TestCredibleSetGate(unittest.TestCase):
    def test_failed_sets_are_dropped_from_both_tables(self):
        # text spellings of the flag, as a TSV read without type inference would give
        arch, members = fixture(
            gate=["TRUE", "FALSE", "True"],
            masses=[(0.2, 0.0, 0.1, 0.7), (np.nan,) * 4, (0.0, 0.0, 0.0, 1.0)],
            member_pips=[[0.5, 0.5], [np.nan], [1.0]])
        kept_arch, kept_members, excluded = m.apply_gate(arch, members)
        self.assertEqual(sorted(kept_arch.credible_set_uid), ["cs0", "cs2"])
        self.assertEqual(sorted(kept_members.credible_set_uid.unique()), ["cs0", "cs2"])
        self.assertEqual(excluded.set_index("credible_set_uid")["reason"].to_dict(),
                         {"cs1": "pip_sum_gate_failed"})

    def test_nonfinite_mass_is_excluded_not_summed_as_zero(self):
        arch, members = fixture(
            gate=[True, True, True, True],
            masses=[(0.2, 0.0, 0.1, 0.7), (np.nan, 0.0, 0.5, 0.5),
                    (0.0, 0.0, 0.0, 1.0), (0.0, 0.0, 0.0, 1.0)],
            member_pips=[[1.0], [1.0], [0.4, np.inf], [1.0]])
        kept_arch, kept_members, excluded = m.apply_gate(arch, members)
        self.assertEqual(sorted(kept_arch.credible_set_uid), ["cs0", "cs3"])
        self.assertEqual(sorted(kept_members.credible_set_uid.unique()), ["cs0", "cs3"])
        self.assertEqual(excluded.set_index("credible_set_uid")["reason"].to_dict(),
                         {"cs1": "nonfinite_architecture_mass",
                          "cs2": "nonfinite_member_pip"})

    def test_disagreeing_gate_flags_raise(self):
        arch, members = fixture(gate=[True, True], masses=[(0, 0, 0, 1)] * 2,
                                member_pips=[[1.0], [1.0]])
        members.loc[members.credible_set_uid == "cs1", "pip_sum_gate_passed"] = False
        with self.assertRaises(ValueError):
            m.apply_gate(arch, members)


if __name__ == "__main__":
    unittest.main()
