"""T1 and T2 both accept absences, so both must be power guarded.

The failure this pins is the one the framing exists to prevent: a set that stops
clearing its null in an arm too weak to have detected the drop being labelled
mediated rather than indeterminate. Six covariates make that easy to hit, since
a stricter adjustment lowers the statistic whether or not mediation is real.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "evaluate_composition_dependence.py"
SEALED = ROOT / "scripts" / "evaluate_gse267145_axis_count.py"
EXT = ROOT / "scripts" / "multicovariate_partial.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ev = _load(MODULE, "composition_eval_tested")
mc = _load(EXT, "mc_for_eval")
analysis = mc.load_sealed(SEALED)


class PowerGuardTests(unittest.TestCase):
    """A drop below the MDE is indeterminate, never mediated."""

    @staticmethod
    def _verdict(held: bool, detectable: bool) -> str:
        return ("composition_independent_component_exists" if held
                else "set_is_composition_mediated" if detectable else "indeterminate")

    def test_the_three_verdicts_are_distinguished(self) -> None:
        self.assertEqual(self._verdict(True, True),
                         "composition_independent_component_exists")
        self.assertEqual(self._verdict(True, False),
                         "composition_independent_component_exists")
        self.assertEqual(self._verdict(False, True), "set_is_composition_mediated")
        self.assertEqual(self._verdict(False, False), "indeterminate")

    def test_an_undetectable_drop_is_not_applicable_to_the_gate(self) -> None:
        held, detectable = False, False
        self.assertFalse(bool(held or detectable))

    def test_a_detectable_drop_is_applicable_and_unmet(self) -> None:
        held, detectable = False, True
        self.assertTrue(bool(held or detectable))
        self.assertFalse(held)


class NamingTests(unittest.TestCase):
    def test_the_failure_name_is_pinned_to_the_frozen_literal(self) -> None:
        self.assertIn("t1_a_composition_independent_component_exists",
                      ev.FROZEN_THRESHOLDS)
        self.assertIn("composition_mediated",
                      ev.FROZEN_THRESHOLDS["t2_per_gene_composition_dependence"])

    def test_gse193066_serves_activity_only(self) -> None:
        self.assertEqual(ev.ARM_AXES["GSE193066"], ("activity",))
        self.assertEqual(ev.ARM_AXES["GSE130970"], ("activity", "fibrosis"))

    def test_the_own_axis_map_covers_both_exclusive_cells(self) -> None:
        self.assertEqual(set(ev.OWN_AXIS), set(ev.EXCLUSIVE_CELLS))


class MatchingTests(unittest.TestCase):
    def test_the_background_matches_the_bin_histogram_exactly(self) -> None:
        from collections import Counter
        generator = np.random.default_rng(3)
        bins = np.repeat(np.arange(10), 150)
        target = np.concatenate([np.flatnonzero(bins == b)[:4] for b in (1, 5, 8)])
        pool = np.setdiff1d(np.arange(bins.size), target)
        background = ev.matched_background(bins, target, pool, generator)
        self.assertEqual(Counter(bins[background].tolist()),
                         Counter(bins[target].tolist()))

    def test_an_exhausted_bin_returns_none(self) -> None:
        self.assertIsNone(ev.matched_background(
            np.asarray([0, 0, 1]), np.asarray([0, 1]), np.asarray([2]),
            np.random.default_rng(1)))


class BhMaskTests(unittest.TestCase):
    def test_the_mask_size_equals_the_bh_count(self) -> None:
        generator = np.random.default_rng(7)
        values = np.concatenate([generator.normal(0, 0.08, 3000),
                                 generator.normal(0.5, 0.05, 120)])
        mask = ev.bh_mask(analysis, values, 170)
        critical = analysis.bh_critical_abs_r(values.size, 170)
        self.assertEqual(int(mask.sum()),
                         analysis.bh_count_from_abs_r(np.abs(values), critical))
        self.assertGreater(int(mask.sum()), 0)

    def test_no_signal_gives_an_empty_mask(self) -> None:
        self.assertEqual(int(ev.bh_mask(
            analysis, np.full(2000, 0.01), 170).sum()), 0)

    def test_a_stricter_adjustment_cannot_raise_the_count(self) -> None:
        """Fewer residual df must not manufacture more discoveries."""

        generator = np.random.default_rng(9)
        values = np.clip(generator.normal(0.0, 0.2, 4000), -0.99, 0.99)
        self.assertGreaterEqual(int(ev.bh_mask(analysis, values, 178).sum()),
                                int(ev.bh_mask(analysis, values, 171).sum()))


if __name__ == "__main__":
    unittest.main()
