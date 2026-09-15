"""The build must reproduce Stage 0c exactly and must not touch the outcomes.

Two properties are pinned. The reproduction guard has to be a real abort, not a
warning, because the whole justification for recomputing instead of reading a
digest is that an exact match proves identity. And the four-way assignment must
reconstruct Stage 0c's two published marginals, which is an arithmetic identity
the cell sizes cannot satisfy by accident.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "build_two_axis_external_substrate.py"
ARM_A = ROOT / "scripts" / "evaluate_gse267145_axis_count.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


build = _load(MODULE, "two_axis_substrate_tested")


class InvariantTests(unittest.TestCase):
    def test_the_published_stage_0c_numbers_are_the_pinned_ones(self) -> None:
        inv = build.STAGE_0C_INVARIANTS
        self.assertEqual(inv["realised_universe"], 61940)
        self.assertEqual(inv["nas_score|fibrosis_stage"]["genes_bh_below_0_05"], 4388)
        self.assertEqual(inv["fibrosis_stage|nas_score"]["genes_bh_below_0_05"], 1305)
        self.assertEqual(
            inv["nas_score|fibrosis_stage"]["observed_max_abs_partial_r"], 0.5689)
        self.assertEqual(
            inv["fibrosis_stage|nas_score"]["observed_max_abs_partial_r"], 0.6177)

    def test_the_four_cells_are_exhaustive_and_disjoint(self) -> None:
        self.assertEqual(
            set(build.CELLS),
            {"activity_only", "fibrosis_only", "both", "neither"},
        )
        self.assertEqual(len(build.CELLS), 4)

    def test_the_marginal_identity_the_build_asserts(self) -> None:
        """activity_only + both must equal 4388; fibrosis_only + both must equal 1305."""

        activity = np.zeros(100, dtype=bool)
        fibrosis = np.zeros(100, dtype=bool)
        activity[:40] = True
        fibrosis[30:45] = True
        assignment = np.where(
            activity & fibrosis, "both",
            np.where(activity, "activity_only",
                     np.where(fibrosis, "fibrosis_only", "neither")))
        sizes = {c: int((assignment == c).sum()) for c in build.CELLS}
        self.assertEqual(sizes["activity_only"] + sizes["both"], int(activity.sum()))
        self.assertEqual(sizes["fibrosis_only"] + sizes["both"], int(fibrosis.sum()))
        self.assertEqual(sum(sizes.values()), 100)


class NormalizationTests(unittest.TestCase):
    def test_cpm_makes_every_sample_sum_to_a_million(self) -> None:
        values = np.asarray([[1.0, 2.0], [3.0, 8.0], [6.0, 10.0]])
        cpm = build.counts_per_million(values)
        for total in cpm.sum(axis=0):
            self.assertAlmostEqual(float(total), 1e6, places=6)

    def test_a_zero_library_is_refused_rather_than_dividing(self) -> None:
        with self.assertRaises(build.SubstrateError):
            build.counts_per_million(np.asarray([[0.0, 1.0], [0.0, 2.0]]))

    def test_cpm_preserves_within_sample_gene_order(self) -> None:
        values = np.asarray([[5.0], [1.0], [9.0]])
        cpm = build.counts_per_million(values)
        self.assertEqual(list(np.argsort(values[:, 0])), list(np.argsort(cpm[:, 0])))


class MembershipTests(unittest.TestCase):
    def test_membership_size_equals_the_bh_count(self) -> None:
        analysis = _load(ARM_A, "arm_a_for_membership")
        generator = np.random.default_rng(4)
        n, g = 60, 900
        exposure = generator.integers(0, 8, size=n).astype(float)
        block = generator.normal(size=(n, g))
        block[:, :80] += 1.4 * (exposure - exposure.mean())[:, None]
        ranks = analysis.average_ranks(block)
        ranks -= ranks.mean(axis=0, keepdims=True)
        covariate = analysis.unit_ranks(generator.integers(0, 5, size=n).astype(float))
        scaled, _ = analysis.unit_columns(analysis.residualize(ranks, covariate))
        member, correlations, count = build.bh_membership(
            analysis, scaled,
            analysis.residualize(analysis.centred_ranks(exposure), covariate), n)
        self.assertEqual(int(member.sum()), count)
        self.assertGreater(count, 0)

    def test_no_signal_gives_an_empty_membership_not_a_full_one(self) -> None:
        analysis = _load(ARM_A, "arm_a_for_empty_membership")
        generator = np.random.default_rng(5)
        n, g = 60, 900
        ranks = analysis.average_ranks(generator.normal(size=(n, g)))
        ranks -= ranks.mean(axis=0, keepdims=True)
        covariate = analysis.unit_ranks(generator.integers(0, 5, size=n).astype(float))
        scaled, _ = analysis.unit_columns(analysis.residualize(ranks, covariate))
        exposure = generator.integers(0, 8, size=n).astype(float)
        member, _, count = build.bh_membership(
            analysis, scaled,
            analysis.residualize(analysis.centred_ranks(exposure), covariate), n)
        self.assertEqual(count, int(member.sum()))
        self.assertLess(count, g // 2)


if __name__ == "__main__":
    unittest.main()
